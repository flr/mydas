require(data.table)


# readFLSss3 {{{

#' A function to read SS3 results as an FLStock object
#'
#' Results of a run of the Stock Synthesis sofware, SS3 (Methot & Wetzel, 2013),
#' can be loaded into an object of class \code{\link{FLStock}}. The code makes
#' use of the r4ss::SS_output function to obtain a list from Report.sso. The
#' following elements of that list are used to generate the necessary information
#' for the slots in \code{\link{FLStock}}: "catage", "natage", "ageselex",
#' "endgrowth", "catch_units", "nsexes", "nseasons", "nareas", "IsFishFleet",
#' "fleet_ID", "FleetNames", "spawnseas", "inputs" and "SS_version".
#'
#' @references
#' Methot RD Jr, Wetzel CR (2013) Stock Synthesis: A biological and statistical
#' framework for fish stock assessment and fishery management.
#' Fisheries Research 142: 86-99.
#'
#' @param dir Directory holding the SS3 output files
#' @param birthseas Birth seasons for this stock, defaults to spawnseas
#' @param name Name of the output object to fil the name slot
#' @param desc Description of the output object to fill the desc slot
#' @param ... Any other argument to be passed to `r4ss::SS_output`
#'
#' @return An object of class `\link{FLStock}`
#'
#' @name readFLSss3
#' @rdname readFLSss3
#' @aliases readFLSss3
#'
#' @author The FLR Team
#' @seealso \link{FLComp}
#' @keywords classes

readFLSss3 <- function(dir, repfile="Report.sso", compfile="CompReport.sso",
                       ...) {
  
  out <- readOutputss3(dir, repfile=repfile, compfile=compfile)
  
  if(out$SS_versionNumeric > 3.24)
    res <- buildFLSss330(out, ...)
  else
    res <- buildFLSss3(out, ...)
  
  # CHANGE mat and *.wt if wtatage file being used 
  if(out$wtatage_switch) {
    
    # FIND wtatage.ss_new
    waafile <- list.files(dir)[grep("wtatage.ss_new", list.files(dir))]
    
    # LOAD wtatage.ss_new
    waa <- data.table(SS_readwtatage(file.path(dir, waafile)))
    
    # SET year, unit and season
    waa[, year:=abs(Yr)]
    waa[, unit:=Sex]
    waa[, season:=Seas]
    
    # GET ages
    ages <- dimnames(res)$age
    
    # SUBSET FLStock years
    waa <- waa[year %in% dimnames(res)$year,]
    
    # SPLIT weights by fleet
    was <- split(waa, by="Fleet")
    
    # CREATE FLQuants
    wasq <- lapply(was, function(x)
      as.FLQuant(melt(x[, -seq(1, 6)], id=c("unit", "year", "season"),
                      measure=ages, variable.name = "age", value.name = "data")))
    
    # stock.wt, Fleet = 0
    stock.wt(res)[] <- wasq[["0"]]
    
    # mat, Fleet = -2 / wt
    nmat <- wasq[["-2"]] / wasq[["0"]]
    mat(res)[] <- nmat
    
    # IDENTIFY catch fleets
    if(is.null(out$fleet_type)) {
      out$fleet_type <- rep(3, out$nfleets)
      out$fleet_type[out$fleet_ID %in% unique(out$catch$Fleet)] <- 1
    }
    
    idx <- names(wasq)[!names(wasq) %in% c("0", "-1", "-2")][out$fleet_type == 1]
    
    # COMPUTE catch.wt DEBUG weighted average
    catch.wt(res)[] <- Reduce("+", wasq[idx]) /
      (length(idx))
    landings.wt(res) <- catch.wt(res)
    discards.wt(res) <- catch.wt(res)
    
    catch(res) <- computeCatch(res)
    landings(res) <- computeLandings(res)
    discards(res) <- computeDiscards(res)
    stock(res) <- computeStock(res)
  }
  
  return(res)
  
} # }}}


buildFLSss330 <- function(out, birthseas=out$birthseas, name=out$Control_File,
                          desc=paste(out$inputs$repfile, out$SS_versionshort, sep=" - "),
                          fleets=setNames(nm=out$fleet_ID[out$IsFishFleet]), range="missing") {
  
  # DIMENSIONS
  dims <- dimss3(out)
  
  # SUBSET out
  out <- out[c("catage", "natage", "ageselex", "endgrowth", "Control_File",
               "catch_units", "nsexes", "nseasons", "nareas", "IsFishFleet", "fleet_ID",
               "FleetNames", "birthseas", "spawnseas", "inputs", "SS_versionshort",
               "discard", "discard_at_age", "catch", "NatMort_option", "GrowthModel_option",
               "Maturity_option", "Fecundity_option", "Z_at_age", "M_at_age",
               "mean_body_wt", "Spawn_timing_in_season")]
  
  # GET ages from catage
  ages <- getRange(out$catage)
  ages <- ac(seq(ages['min'], ages['max']))
  dmns <- getDimnames(out, birthseas=birthseas)
  dim <- unlist(lapply(dmns, length))
  
  # ENDGROWTH
  if(out$nsexes == 1) {
    endgrowth <- data.table(out$endgrowth,
                            key=c("Seas", "Platoon", "Settlement", "int_Age"))
  } else {
    endgrowth <- data.table(out$endgrowth,
                            key=c("Seas", "Sex", "Platoon", "Settlement", "int_Age"))
  }
  
  # SET Age and unit
  endgrowth[, Age:=int_Age]
  endgrowth[, unit:=codeUnit(Sex, Platoon)]
  
  # NATAGE
  natage <- data.table(out$natage)
  natage[, unit:=codeUnit(Sex, Platoon)]
  
  # CATCH.N
  catage <- data.table(out$catage)
  catage[, unit:=codeUnit(Sex)]
  # NOTE catage$0 comes out as integer
  catage[, `0` := as.double(`0`)]
  setkey(catage, "Area", "Fleet", "unit", "Yr", "Seas", "Era")
  
  # WT
  wtatage <- endgrowth[, c("Seas", "unit", "Age", paste0("RetWt:_", fleets)),
                       with=FALSE]
  
  # STOCK.WT
  wt <- ss3wt30(endgrowth, dmns, birthseas=1)
  
  # MAT
  mat <- ss3mat30(endgrowth, dmns, birthseas, option=out$Maturity_option)
  
  # CORRECT Mat*Fecund to by unit body weight
  if(out$Maturity_option == 6)
    mat <- mat / wt
  
  # M
  m <- ss3m30(endgrowth, dmns, birthseas)
  
  # STOCK.N
  n <- ss3n30(natage, dmns, birthseas)
  
  # CATCH 
  catches <- ss3catch30(catage, wtatage, dmns, birthseas, fleets)
  
  # CALCULATE total catch.n
  catch.n <- FLQuant(0, dimnames=dmns, units="1000")
  
  for (i in seq(length(fleets)))
    catch.n <- catch.n %++% catches[[i]]$catch.n
  
  # AVERAGE catch.wt weighted by catch.n
  catch.wt <-  FLCore::expand(Reduce("+", lapply(catches,
                                                 function(x) x$catch.n %*% x$catch.wt)) %/% catch.n,
                              year=dmns$year, area=dmns$area)
  
  catch.wt[is.na(catch.wt)] <- (Reduce("+",
                                       lapply(catches, '[[', 'catch.wt')) / length(catches))[is.na(catch.wt)]
  
  # DISCARDS
  if(!is.na(out["discard"])) {
    
    stop("PROBLEM with discards")
    
    datage <- data.table(out$discard_at_age)
    datage[, unit:=codeUnit(Sex)]
    setkey(datage, "Area", "Fleet", "unit", "Yr", "Seas", "Era", "Type")
    
    # DEBUG
    ageselex <- data.table(out$ageselex)
    lastyr <- unique(ageselex[Factor=="Asel2", Yr])
    
    # TOTAL selex (catch$kill_nums)
    seltot <- ss3sel.pattern(ageselex, lastyr, fleets, morphs=unique(ageselex$Morph),
                             factor="sel_nums")
    
    # RETAINED selex (catch$ret_num)
    selret <- ss3sel.pattern(ageselex, lastyr, fleets, morphs=unique(ageselex$Morph),
                             factor="sel*ret_nums")
    
    # DEAD selex (catch$kill_num)
    seldea <- ss3sel.pattern(ageselex, lastyr, fleets, morphs=unique(ageselex$Morph),
                             factor="dead_nums")
    
    # DISCARD selex (dead + alive)
    seldis <- mapply(function(x, y) x - y, seltot, selret, SIMPLIFY=FALSE)
    
    discard <- data.table(out$discard)
    catch <- data.table(out$catch)
    
    # F_discards by fleet: catch from last estimation yr
    Fdiscards <- FLQuants(mapply(function(x, y) x * y, seldis,
                                 as.list(catch[Yr == max(Yr) & Fleet %in% fleets, F]), SIMPLIFY=FALSE))
    
    # APPLY Baranov for discards.n
    discards.n <- Reduce('+', mapply(function(x)
      FLQuant((x %/% (x %+% m)) * (1 - exp(-(x %+% m))) * n, units="1000"),
      Fdiscards, SIMPLIFY=FALSE))
    dimnames(discards.n) <- dimnames(catch.n)
    
    # SET discards.n not in discard period as 0
    discards.n[, !dimnames(discards.n)$year %in% discard$Yr] <- 0
    
    discards.wt <- catch.wt
    
  } else {
    discards.n <- catch.n * 0
    discards.wt <- catch.wt
  }
  
  # FLStock
  stock <- FLStock(
    name=name, desc=desc,
    catch.n=catch.n, catch.wt=catch.wt,
    discards.n=discards.n, discards.wt=discards.wt,
    landings.n=catch.n - discards.n, landings.wt=catch.wt,
    stock.n=n, stock.wt=wt,
    m=m, mat=mat)
  
  # CALCULATE stock, catch, landings & discards
  landings(stock) <- computeLandings(stock)
  discards(stock) <- computeDiscards(stock)
  catch(stock) <- computeCatch(stock, slot='all')
  stock(stock) <- computeStock(stock)
  
  # ASSIGN harvest.spwn and m.spwn
  harvest.spwn(stock) <- out$Spawn_timing_in_season
  m.spwn(stock) <- out$Spawn_timing_in_season
  
  # HARVEST
  harvest(stock) <- harvest(stock.n(stock), catch=catch.n(stock), m=m(stock))
  
  # range
  if(!missing(range))
    range(stock) <- range
  
  return(stock)
  
} # }}}



buildFLSss3 <- function(out, birthseas=out$birthseas, name=out$Control_File,
                        desc=paste(out$inputs$repfile, out$SS_versionshort, sep=" - "), range="missing",
                        fleets=setNames(out$fleet_ID[out$IsFishFleet], out$fleet_ID[out$IsFishFleet])) {
  
  # SUBSET out
  out <- out[c("catage", "natage", "ageselex", "endgrowth", "Control_File",
               "catch_units", "nsexes", "nseasons", "nareas", "IsFishFleet", "fleet_ID",
               "FleetNames", "birthseas", "spawnseas", "inputs", "SS_versionshort",
               "discard", "catch")]
  
  # TODO: call spread()
  
  # GET ages from catage
  ages <- getRange(out$catage)
  ages <- ac(seq(ages['min'], ages['max']))
  dmns <- getDimnames(out, birthseas=birthseas)
  dim <- unlist(lapply(dmns, length))
  
  # EXTRACT from out
  if(out$nsexes == 1) {
    endgrowth <- data.table(out$endgrowth, key=c("Seas", "Age"))
  } else {
    endgrowth <- data.table(out$endgrowth, key=c("Seas", "Sex", "Age"))
  }
  # NATAGE
  natage <- data.table(out$natage)
  
  # CATCH.N
  catage <- data.table(out$catage)
  setkey(catage, "Area", "Fleet", "Gender", "Morph", "Yr", "Seas", "Era")
  
  # STOCK.WT
  wt <- ss3wt(endgrowth, dmns, birthseas)
  
  # MAT
  mat <- ss3mat(endgrowth, dmns, birthseas)
  
  # M
  m <- ss3m(endgrowth, dmns, birthseas)
  
  # STOCK.N
  n <- ss3n(natage, dmns, birthseas)
  
  # CATCH.WT, assumes _mat_option == 3
  wtatage <- endgrowth[BirthSeas %in% birthseas,
                       c("Seas", "Sex", "BirthSeas", "Age", paste0("RetWt:_", fleets)), with=FALSE]
  
  catches <- ss3catch(catage, wtatage, dmns, birthseas, fleets)
  
  # CALCULATE total catch.n & catch.wt
  catch.n <- FLQuant(0, dimnames=dmns, units="1000")
  catch.wt <- FLQuant(0, dimnames=dmns, units="kg")
  catch.den <- FLQuant(0, dimnames=dmns, units="kg")
  
  for (i in seq(length(fleets))) {
    catch.n <- catch.n %++% catches[[i]]$catch.n
    
    # DROP NAs from discards fleets
    if(!all(is.na(catches[[i]]$catch.wt))) {
      catch.den <- catch.n %++% catches[[i]]$catch.n
      catch.wt <- catch.wt %++% (catches[[i]]$catch.wt * (catches[[i]]$catch.n + 1e-3))
    }
  }
  
  # COMPUTE no. of fleets per area
  wtsbyarea <- c(table(unlist(lapply(catches, function(x) dimnames(x$catch.wt)$area))))
  
  # DIVIDE by total catch
  catch.wt <- catch.wt / (catch.den + FLQuant(rep(1e-3 * wtsbyarea,
                                                  each=prod(dim(catch.n)[-5])), dimnames=dmns, units="1000"))
  
  # RESET units(catch.wt)
  units(catch.wt) <- "kg"
  
  # DISCARDS
  if(!is.na(out["discard"])) {
    
    ageselex <- data.table(out$ageselex)
    lastyr <- unique(ageselex[Factor=="sel_nums", Yr])
    lastyr <- unique(ageselex[Factor=="Asel2", Yr])
    
    # TOTAL selex (catch$kill_nums)
    seltot <- ss3sel.pattern(ageselex, lastyr, fleets, morphs=unique(ageselex$Morph),
                             factor="sel_nums")
    
    # RETAINED selex (catch$ret_num)
    selret <- ss3sel.pattern(ageselex, lastyr, fleets, morphs=unique(ageselex$Morph),
                             factor="sel*ret_nums")
    
    # DEAD selex (catch$kill_num)
    seldea <- ss3sel.pattern(ageselex, lastyr, fleets, morphs=unique(ageselex$Morph),
                             factor="dead_nums")
    
    # DISCARD selex (dead + alive)
    seldis <- mapply(function(x, y) x - y, seltot, selret, SIMPLIFY=FALSE)
    
    discard <- data.table(out$discard)
    catch <- data.table(out$catch)
    
    # F_discards by fleet: catch from last estimation yr
    Fdiscards <- FLQuants(mapply(function(x, y) x * y, seldis,
                                 as.list(catch[Yr == max(Yr) & Fleet %in% fleets, F]), SIMPLIFY=FALSE))
    
    # APPLY Baranov for discards.n
    discards.n <- Reduce('+', mapply(function(x)
      FLQuant((x %/% (x %+% m)) * (1 - exp(-(x %+% m))) * n, units="1000"),
      Fdiscards, SIMPLIFY=FALSE))
    dimnames(discards.n) <- dimnames(catch.n)
    
    # SET discards.n not in discard period as 0
    discards.n[, !dimnames(discards.n)$year %in% discard$Yr] <- 0
    
    discards.wt <- catch.wt
    
  } else {
    discards.n <- catch.n
    discards.n[] <- 0
    discards.wt <- catch.wt
  }
  
  # FLStock
  stock <- FLStock(
    name=name, desc=desc,
    catch.n=catch.n, catch.wt=catch.wt,
    discards.n=discards.n, discards.wt=discards.wt,
    landings.n=catch.n - discards.n, landings.wt=catch.wt,
    stock.n=n, stock.wt=wt,
    m=m, mat=mat)
  
  # CALCULATE stock, catch, landings & discards
  landings(stock) <- computeLandings(stock)
  discards(stock) <- computeDiscards(stock)
  catch(stock) <- computeCatch(stock, slot='all')
  stock(stock) <- computeStock(stock)
  
  # ASSIGN harvest.spwn and m.spwn in birthseas
  harvest.spwn(stock)[,,,birthseas] <- 0
  m.spwn(stock)[,,,birthseas] <- 0
  
  # HARVEST
  harvest(stock) <- harvest(stock.n(stock), catch=catch.n(stock), m=m(stock))
  
  # range
  if(!missing(range))
    range(stock) <- range
  
  return(stock)
  
} # }}}

readOutputss3 <- function(dir, repfile = "Report.sso",
                          compfile = "CompReport.sso", compress="gz") {
  
  # Possibly compressed files
  cfiles <- c(repfile = repfile, compfile = compfile)
  
  # CHECK compressed files
  idx <- file.exists(file.path(dir, paste(cfiles, compress, sep = ".")))
  cfiles[idx] <- paste(cfiles, compress, sep = ".")
  
  out <- SS_output(dir, verbose=FALSE, hidewarn=TRUE, warn=FALSE,
                   printstats=FALSE, covar=FALSE, forecast=FALSE,
                   repfile=cfiles["repfile"], compfile=cfiles["compfile"])
  
  return(out) }

dimss3 <- function(out) {
  
  range <- getRange(out$catage)
  
  return(list(
    age = length(seq(range["min"], range["max"])),
    year   = length(seq(range["minyear"], range["maxyear"])),
    sex = out$nsexes,
    biop = out$N_bio_patterns,
    gpatterns = out$ngpatterns,
    platoon = out$N_platoons,
    # unit = sex * platoon
    unit = out$nsexes * out$N_platoons,
    season = out$nseasons,
    spwn = (out$Spawn_seas - 1) * (1 / out$nseasons) +
      out$Spawn_timing_in_season * (1 / out$nseasons),
    subseason = out$N_sub_seasons,
    area = out$nareas,
    fleet = sum(out$IsFishFleet),
    index = sum(!out$IsFishFleet),
    M_option = out$NatMort_option,
    growth_option = out$GrowthModel_option,
    mat_option = out$Maturity_option,
    fec_option = out$Fecundity_option
  ))}

getRange <- function(x) {
  
  # empty range
  range <- rep(as.numeric(NA), 7)
  names(range) <- c("min", "max", "plusgroup", "minyear", "maxyear",
                    "minfbar", "maxfbar")
  
  # age range from catage
  # TODO FIND more secure way to find ages columns
  idx <- grep("Era", names(x))
  range[c("min", "max")] <- range(as.numeric(names(x)[-seq(1, idx)]))
  
  # plusgroup = max
  range["plusgroup"] <- range["max"]
  
  # min/maxfbar = min/max
  range[c("minfbar", "maxfbar")] <- range[c("min", "max")]
  
  # year range from catage
  range[c("minyear", "maxyear")] <- range(x$Yr[x$Era == "TIME"])
  
  return(range)} 

getDimnames <- function(out, birthseas) {
  
  # GET range
  range <- getRange(out$catage)
  ages <- ac(seq(range['min'], range['max']))
  
  dmns <- list(age=ages,
               year=seq(range['minyear'], range['maxyear']),
               # unit = combinations(Sex, birthseas)
               unit=c(t(outer(switch(out$nsexes, "", c("F", "M")),
                              switch((length(birthseas) > 1) + 1, "", birthseas), paste0))),
               season=switch(ac(out$nseasons), "1"="all", seq(out$nseasons)),
               area=switch(ac(out$nareas), "1"="unique", seq(out$nareas)),
               iter=1)
  
  # TODO HACK
  if(all(dmns$unit == ""))
    dmns$unit <- "unique"
  
  return(dmns)}

codeUnit <- function(Sex, Platoon="missing") {
  
  # SEX as "F" / "M"
  Sex <- if(length(unique(Sex)) == 1){""} else {c("F","M")[Sex]}
  
  # Platoon
  if(!missing(Platoon))
    Platoon <- if(length(unique(Platoon)) == 1){""} else {Platoon}
  else
    Platoon <- ""
  
  unit <- ifelse(paste0(Sex, Platoon) == "", "unique", paste0(Sex, Platoon))
  
  return(unit)}

# ss3slots.R - DESC
# ss3om/R/ss3slots.R

# Copyright European Union, 2015-2019; WMR, 2020.
# Author: Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

#' @rdname ss3slot
#' @aliases ss3mat30
#' @details - `ss3mat30` returns the `mat` slot.

ss3mat30 <- function(endgrowth, dmns, birthseas, option=3) {
  
  # EXTRACT mat - endgrowth
  mat <- endgrowth[, .(unit, Seas, Age, Age_Mat, `Mat*Fecund`, Wt_Beg,
                       Mat_F_wtatage, Mat_F_Natage)]
  
  # maturity option 6: mat=Mat*Fecund / max(Mat*Fecund)
  if(option == 6)
    #mat[, mat:= `Mat*Fecund` / max(`Mat*Fecund`), by=.(unit, Seas)]
    mat[, mat:= `Mat*Fecund`]
  
  # maturity option 3: mat=Age_Mat
  if(option == 3)
    mat[, mat:= Age_Mat]
  
  # maturity option 1: mat=Mat*Fecund / Wt_Beg
  if(option == 1)
    mat[, mat:= `Mat*Fecund` / Wt_Beg]
  
  # DEBUG
  if(option == 2)
    mat[, mat:= Age_Mat]
  
  # DEBUG
  if(option == 5)
    mat[, mat:= 0]
  
  # DEBUG
  if(!option %in% c(3, 6, 1, 2, 5))
    stop(paste0("maturity option not covered yet, option: ", option))
  
  # DELETE columns
  mat[ ,`:=`(Age_Mat = NULL, `Mat*Fecund` = NULL, Wt_Beg = NULL,
             Mat_F_wtatage = NULL, Mat_F_Natage = NULL)]
  
  # RENAME
  names(mat) <- c("unit", "season", "age", "data")
  
  # TURN -1/NaN to 0
  mat[, data:=ifelse(is.nan(data), 0, data)]
  mat[, data:=ifelse(data==-1, 0, data)]
  
  # EXPAND by year & unit & area
  mat <- FLCore::expand(as.FLQuant(mat[, .(season, unit, age, data)],
                                   units=""), year=dmns$year, unit=dmns$unit, season=dmns$season, area=dmns$area)
  
  return(mat)
}

#' @rdname ss3slot
#' @aliases ss3m
#' @details - `ss3m` returns the `m` slot.

ss3m30 <- function(endgrowth, dmns, birthseas) {
  
  # EXTRACT m - biol[, Seas, BirthSeas, Age, M]
  m <- endgrowth[, .(Age, unit, Seas, M)]
  
  # TODO SPLIT M across seasons
  m[, M:=M/length(dmns$season)]
  
  # RENAME
  names(m) <- c("age", "unit", "season", "data")
  
  # EXPAND by year, unit, season & area
  # BUG expand not filling
  m <- FLCore::expand(as.FLQuant(m[,.(season, age, data, unit)], units="m"),
                      year=dmns$year, unit=dmns$unit, season=dmns$season, area=dmns$area)
  
  return(m)
}

#' @rdname ss3slot
#' @aliases ss3n
#' @param n A data frame obtained from SS_output$natage.
#' @details - `ss3n30` returns the `stock.n` slot.

ss3n30 <- function(n, dmns, birthseas) {
  
  # SELECT start of season (Beg/Mid == 'B'), Era == 'TIME' & cols
  n <- n[`Beg/Mid` == "B" & Era == 'TIME',
         .SD, .SDcols = c("Area", "unit", "Yr", "Seas", dmns$age)]
  
  # MELT by Sex, BirthSeas, Yr & Seas
  n <- data.table::melt(n, id.vars=c("Area", "unit", "Yr", "Seas"),
                        variable.name="age")
  
  # SUBSET according to birthseas
  # n <- n[BirthSeas %in% birthseas,]
  
  # RENAME
  names(n) <- c("area", "unit", "year", "season", "age", "data")
  n <- as.FLQuant(n, units="1000")
  dimnames(n) <- dmns
  
  return(n)
}

#' @rdname ss3slot
#' @aliases ss3catch
#' @param catage A data frame obtained from SS_output$catage.
#' @param wtatage A data frame obtained from SS_output$endgrowth but subset for `birthseas` and `RetWt:_idx`.
#' @param idx The fishing fleets, as in `SS_output$fleet_ID[SS_output$IsFishFleet]`.
#' @details - `ss3catch` currently returns the `landings.n` slot, equal to `catch.n` as discards are not being parsed.

ss3catch30 <- function(catage, wtatage, dmns, birthseas, idx) {
  
  # FIND and SUBSET fishing fleets, TIME and BirthSeas
  catage <- catage[Fleet %in% idx & Era == "TIME",]
  
  # RENAME Area and Season if only 1
  cols <- c("Seas", "Area")
  catage[, (cols) := lapply(.SD, as.character), .SDcols = cols]
  catage[, Seas := if(length(unique(Seas)) == 1) "all" else Seas]
  catage[, Area := if(length(unique(Area)) == 1) "unique" else Area]
  
  # MELT by Sex, BirthSeas, Yr & Seas
  catage <- data.table::melt(catage, id.vars=c("Area", "Fleet", "Yr", "Seas", "unit"),
                             measure.vars=dmns$age, variable.name="age")
  names(catage) <- c("area", "fleet", "year", "season", "unit", "age", "data")
  
  # RENAME Area and Season if only 1
  cols <- c("Seas")
  wtatage[, (cols) := lapply(.SD, as.character), .SDcols = cols]
  wtatage[, Seas := if(length(unique(Seas)) == 1) "all" else Seas]
  
  # MELT by Sex, BirthSeas, Yr & Seas
  wtatage <- data.table::melt(wtatage, id.vars=c("Age", "unit", "Seas"),
                              measure.vars=paste0("RetWt:_", idx), variable.name="fleet")
  names(wtatage) <- c("age", "unit", "season", "fleet", "data")
  wtatage[,fleet:=sub("RetWt:_", "", fleet)]
  
  # FLQuants for catch per fleet
  catch <- lapply(idx, function(x) {
    catch.n <- as.FLQuant(catage[fleet %in% x,][, fleet:=NULL], units="1000")
    catch.wt <- do.call('expand',
                        c(list(x=as.FLQuant(wtatage[fleet %in% x,][, fleet:=NULL], units="kg")),
                          dimnames(catch.n)[c("year", "area")]))
    return(FLQuants(catch.n=catch.n, catch.wt=catch.wt))
  }
  )
  return(catch)
} 

ss3z30 <- function(zatage, m, dmns) {
  
  zaa <- zatage[Yr %in% dmns$year, -1]
  setnames(zaa, c("Sex", "Yr"), c("unit", "year"))
  
  zatage <- data.table::melt(zaa, id.vars=c("unit", "year"),
                             measure.vars=dmns$age, variable.name="age", value.name = "data")
  
  z <- as.FLQuant(zatage, units="z")
  dimnames(z) <- dmns[-4]
  
  return(z)
}

ss3wt30 <- function(endgrowth, dmns, birthseas) {
  
  # EXTRACT
  wt <- endgrowth[, list(Age, unit, Seas, Wt_Beg)]
  
  # RENAME
  names(wt) <- c("age", "unit", "season", "data")
  
  # EXPAND by year, unit & season
  return(FLCore::expand(as.FLQuant(wt, units="kg"),
                        year=dmns$year, unit=dmns$unit, season=dmns$season, area=dmns$area))
}
