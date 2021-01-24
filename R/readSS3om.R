library(data.table)

ss3wt <- function(endgrowth, dmns, birthseas) {
  
  # EXTRACT stock.wt - Wt_mid @ endgrowth[, Seas, BirthSeas, Age, M]
  wt <- endgrowth[BirthSeas %in% birthseas,
                  list(BirthSeas, Sex, Seas, Age, Wt_Beg)]
  
  # CREATE unit from Sex + BirthSeas
  wt[, uSex:={if(length(unique(Sex)) == 1){""} else {c("F","M")[Sex]}}]
  wt[, uBirthSeas:={if(length(unique(BirthSeas)) == 1){""} else {BirthSeas}}]
  wt[, unit:=paste0(uSex, uBirthSeas)]
  wt[, unit:=ifelse(paste0(uSex, uBirthSeas) == "", "unique", paste0(uSex, uBirthSeas))]
  wt[, c("Sex","uSex","BirthSeas","uBirthSeas") := NULL]
  
  # RENAME
  names(wt) <- c("season", "age", "data", "unit")
  
  # EXPAND by year, unit & season
  return(FLCore::expand(as.FLQuant(wt[, .(season, age, data, unit)], units="kg"),
                        year=dmns$year, unit=dmns$unit, season=dmns$season, area=dmns$area))}

ss3mat <- function(endgrowth, dmns, birthseas, option=3) {
  
  # EXTRACT mat - endgrowth
  mat <- endgrowth[BirthSeas %in% birthseas,
                   # Mat_Numbers
                   list(BirthSeas, Sex, Seas, Age, Age_Mat, `Mat*Fecund`, Wt_Beg)]
  
  # TODO CHECK maturity_option 2, 4, 5
  # IF maturity_option == 3, mat = mat / wt
  if(all(mat[, Age_Mat] %in% c(0,1)))
    mat[, `Mat*Fecund`:= `Mat*Fecund` / Wt_Beg]
  
  # maturity option 3: mat=Age_Mat
  if(option == 3)
    mat[, `Mat*Fecund`:= Age_Mat]
  
  mat[ ,`:=`(Age_Mat = NULL, Wt_Beg = NULL)]
  
  # RENAME
  names(mat) <- c("BirthSeas", "Sex", "season", "age", "data")
  
  # TURN -1 to 0
  # TODO CHECK M values == 0
  # TODO TURN !birthseas to 0
  mat[, data:=ifelse(data==-1, 0, data)]
  
  # SWT unit from Sex and BirthSeas
  mat[, uSex:={if(length(unique(Sex)) == 1){""} else {c("F","M")[Sex]}}]
  mat[, uBirthSeas:={if(length(unique(BirthSeas)) == 1){""} else {BirthSeas}}]
  mat[, unit:=ifelse(paste0(uSex, uBirthSeas) == "", "unique", paste0(uSex, uBirthSeas))]
  mat[ ,c("Sex","uSex","BirthSeas","uBirthSeas") := NULL]
  
  # EXPAND by year & unit & area
  mat <- FLCore::expand(as.FLQuant(mat[, .(season, unit, age, data)],
                                   units=""), year=dmns$year, unit=dmns$unit, season=dmns$season, area=dmns$area)
  
  return(mat)}

ss3m <- function(endgrowth, dmns, birthseas) {
  
  # EXTRACT m - biol[, Seas, BirthSeas, Age, M]
  m <- endgrowth[BirthSeas %in% birthseas,
                 list(BirthSeas, Sex, Seas, Age, M)]
  
  # CREATE unit from Sex + BirthSeas
  m[, uSex:={if(length(unique(Sex)) == 1){""} else {c("F","M")[Sex]}}]
  m[, uBirthSeas:={if(length(unique(BirthSeas)) == 1){""} else {BirthSeas}}]
  m[, unit:=paste0(uSex, uBirthSeas)]
  m[ ,c("Sex","uSex","BirthSeas","uBirthSeas") := NULL]
  
  # SPLIT M across seasons
  m[, M:=M/length(dmns$season)]
  
  # RENAME
  names(m) <- c("season", "age", "data", "unit")
  
  # EXPAND by year, unit, season & area
  # BUG expand not filling
  m <- FLCore::expand(as.FLQuant(m[,.(season, age, data, unit)], units="m"),
                      year=dmns$year, unit=dmns$unit, season=dmns$season, area=dmns$area)
  
  return(m)}

ss3n <- function(n, dmns, birthseas) {
  
  # SELECT start of season (Beg/Mid == 'B'), Era == 'TIME' & cols
  n <- n[`Beg/Mid` == "B" & Era == 'TIME',
         .SD, .SDcols = c("Area", "Sex", "BirthSeas", "Yr", "Seas", dmns$age)]
  
  # MELT by Sex, BirthSeas, Yr & Seas
  n <- data.table::melt(n, id.vars=c("Area", "Sex", "BirthSeas", "Yr", "Seas"),
                        variable.name="age")
  
  # SUBSET according to birthseas
  n <- n[BirthSeas %in% birthseas,]
  
  # DROP Sex and BirthSeas
  
  # CREATE unit from Sex + BirthSeas
  n[, uSex:={if(length(unique(Sex)) == 1){""} else {c("F","M")[Sex]}}]
  n[, uBirthSeas:={if(length(unique(BirthSeas)) == 1){""} else {BirthSeas}}]
  n[, unit:=paste0(uSex, uBirthSeas)]
  n[ ,c("Sex","uSex","BirthSeas","uBirthSeas") := NULL]
  
  
  # RENAME
  names(n) <- c("area", "year", "season", "age", "data", "unit")
  n <- as.FLQuant(n, units="1000")
  dimnames(n) <- dmns
  
  return(n)}

ss3catch <- function(catage, wtatage, dmns, birthseas, idx) {
  
  # RECONSTRUCT BirthSeas from Morph & Sex
  catage[, BirthSeas := Morph - max(Seas) * (Gender - 1)]
  catage <- catage[BirthSeas %in% birthseas,]
  
  # CREATE unit from Sex + BirthSeas
  catage[, uSex:={if(length(unique(Gender)) == 1){""} else {c("F","M")[Gender]}}]
  catage[, uBirthSeas:={if(length(unique(BirthSeas)) == 1){""} else {BirthSeas}}]
  catage[, unit:=ifelse(paste0(uSex, uBirthSeas) == "", "unique", paste0(uSex, uBirthSeas))]
  
  wtatage[, uSex:={if(length(unique(Sex)) == 1){""} else {c("F","M")[Sex]}}]
  wtatage[, uBirthSeas:={if(length(unique(BirthSeas)) == 1){""} else {BirthSeas}}]
  wtatage[, unit:=ifelse(paste0(uSex, uBirthSeas) == "", "unique", paste0(uSex, uBirthSeas))]
  
  # FIND and SUBSET fishing fleets, TIME and BirthSeas
  catage <- catage[Fleet %in% idx & Era == "TIME" & BirthSeas %in% birthseas,]
  
  catage[ ,c("Gender","uSex","BirthSeas","uBirthSeas") := NULL]
  
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
  })
  
  return(catch)} 

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
  
  # set plusgroup to max age
  range["plusgroup"] <- range["maxyear"]
  
  return(range)}

# getDimnames {{{
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

buildFLSss3 <- function(out, birthseas=out$birthseas, name=out$Control_File,
                        desc=paste(out$inputs$repfile, out$SS_versionshort, sep=" - "),
                        fleets=setNames(out$fleet_ID[out$IsFishFleet], out$fleet_ID[out$IsFishFleet])) {
  
  # SUBSET out
  out <- out[c("catage", "natage", "ageselex", "endgrowth", "Control_File",
               "catch_units", "nsexes", "nseasons", "nareas", "IsFishFleet", "fleet_ID",
               "FleetNames", "birthseas", "spawnseas", "inputs", "SS_versionshort",
               "discard", "catch")]
  
  # TODO: call spread()
  
  # GET range from catage
  range <- getRange(out$catage)
  ages <- ac(seq(range['min'], range['max']))
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
  
  return(stock)}

readSS3om<-function(out, birthseas=unique(subset(out$timeseries,!is.na(SpawnBio))[,"Seas"])) {
  
  if(out$SS_versionNumeric > 3.24)
    stop("ss3om currently only supports SS3 <= 3.24")
  
  buildFLSss3(out, birthseas=birthseas)} 


