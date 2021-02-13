utils::globalVariables(c("objFn","setParams<-","setControl<-","fit"))

#' #http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf

#' mse
#' 
#' @title mseMPJabba 
#' 
#' @description 
#' @author Laurence Kell, Sea++
#'  
#' @name mseMPJabba
#' 
#' @aliases mseMPJabba mseMPJabba-method mseMPJabba,FLStock,FLBRP-method
#' 
#' @param om \code{FLStock} object as the operating model
#' @param eq  blah,blah,blah,...
#' @param mp  blah,blah,blah,...
#' @param ftar 1.0  blah,blah,blah,...
#' @param btrig 0.5  blah,blah,blah,...
#' @param fmin 0.05  blah,blah,blah,...
#' @param blim 0.3  blah,blah,blah,...      
#' @param sr_deviances rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.3),
#' @param u_deviances  =rlnorm(dim(om)[6],FLQuant(0,dimnames=list(year=start:end)),0.2),
#' @param interval 3  blah,blah,blah,...
#' @param start \code{numeric}  default is range(om)["maxyear"]-30
#' @param end \code{numeric}  default is range(om)["maxyear"]-interval
#' @param maxF 1.5 
#'
#' @docType methods
#' 
#' @rdname mseMPJabba
#' 
#' @examples
#' \dontrun{
#' data(pl4)
#' }
#' 
#' @export mseMPJabba
#' 

oemFn<-function(om,maxyear=max(dimnames(om)$year),devu=NULL){
  
  ts=ddply(model.frame(FLQuants(om,
                                index  =function(x) stock(x), #catch(x)/fbar(x),
                                catch  =catch,
                                ssb    =ssb, 
                                biomass=stock),drop=TRUE),.(year),with, 
           data.frame(index  =mean(index,na.rm=T), catch=sum(catch,na.rm=T), 
                      biomass=sum(biomass,na.rm=T),ssb=sum(ssb,na.rm=T)))
  ts$catch[is.na(ts$catch)][]=0.001
  ts=transform(ts, index=index/mean(index))
  
  ts=subset(ts,year<=maxyear)
  
  if (is.null(devu)) return(ts)
  
  transform(merge(ts,devu,all.x=T),index=index*exp(residual))[,-6]}


jabba2biodyn<-function(object, phase=c("b0"=-1,"r"=4,"k"=3,"p"=-2,"q"=2,"sigma"=1),
                       min=0.1,max=10){
  
  res=mpb:::biodyn()
  params(res)[]=object$pars[c("r","K","m","psi"),"Median"]-c(0,0,1,0)
  catch(res)   =as.FLQuant(transmute(object$inputseries$catch,year=Yr,data=Total))
  res@stock    =as.FLQuant(data.frame(year=names(object$timeseries[,"mu","B"]),data=object$timeseries[,"mu","B"]))
  res@name     =paste(object$assessment,object$scenario)
  res@desc     ="coerced from JABBA"
  
  indices=list()
  for (i in dimnames(object$inputseries$cpue)[[2]][-1])
    indices[[i]]=as.FLQuant(data.frame(year=object$inputseries$cpue[,1],
                                       data=object$inputseries$cpue[,i]))
  res@indices=FLQuants(indices)
  
  #bug
  setParams(res)=res@indices      
  setControl(res,min=min,max=max)=res@params
  
  for (i in names(phase[names(phase)%in%dimnames(control(res))[[1]]]))
    control(res)[i,"phase"]=phase[i]
  
  if ("q"%in%names(phase))
    control(res)[grep("q",dimnames(control(res))[[1]]),"phase"]=phase["q"]
  
  if ("sigma"%in%names(phase))
    control(res)[grep("s",dimnames(control(res))[[1]]),"phase"]=phase["sigma"]
  
  return(res)}

#http://ices.dk/sites/pub/Publication%20Reports/Advice/2017/2017/12.04.03.01_Reference_points_for_category_1_and_2.pdf

mseMPJabba<-function(om,eq,sa,
                     sr_deviances,u_deviances,
                     ftar=1.0,btrig=000.005,fmin=000.0005,blim=000.003,
                     start=range(om)["maxyear"]-15,end=range(om)["maxyear"],interval=3,
                     maxF=5,bndTac=c(0.8,1.2),
                     path=""){ 
  
  #if ("FLQuant"%in%is(  u_deviances)) u_deviances  =FLQuants(u_deviances)
  #if ("FLQuant"%in%is(sel_deviances)) sel_deviances=FLQuants(sel_deviances)
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(sr_deviances)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Capacity limits maximum F
  maxF=median(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## Loop round years
  cat("==")
  for (iYr in seq(start,range(om)["maxyear"]-interval,interval)){
    cat(iYr,", ",sep="")
    
    #### Observation Error Model, get catch and single CPUE, at the moment
    dat    =oemFn(om,iYr)

    if (dim(u_deviances)[4]>1){
      u=transmute(merge(dat[,c("year","index")],as.data.frame(u_deviances,drop=T),all.x=T),Yr=year,season=season,Index=index*exp(data))}
    else 
      u=transmute(merge(dat[,c("year","index")],as.data.frame(u_deviances,drop=T),all.x=T),Yr=year,season=1,     Index=index*exp(data))
  
    sa$cpue=cast(u,Yr~season,value="Index")[,seq(dim(u_deviances)[4]+1)] 
    sa$catch=transmute(dat,   Yr=year,Total=catch)
    
    #### Management Procedure
    ## Assessment
    sa$se   =sa$cpue
    sa$se[,-1]=0.1
  
    jb=do.call("build_jabba",sa)
    jb=fit_jabba(jb,
                 init.values=TRUE,
                 init.K=sa$K.prior[1],
                 init.r=sa$r.prior[1]*2,
                 init.q=rep(mean(dat$index/dat$biomass),dim(sa$cpue)[2]-1),
                 ni =5500,
                 nt =1,
                 nb =500,
                 nc =2)
    mp=jabba2biodyn(jb)
    mp=window(mp,end=iYr)
    
    ## Update current year catch
    catch(mp)[,ac(rev(iYr-seq(interval+1)))]=aaply(catch(om)[,ac(rev(iYr-seq(interval+1)))],2,sum)
    catch(mp)[,ac(iYr)]=aaply(catch(om)[,ac(iYr)],2,sum)
    mp=fwd(mp,catch=catch(mp)[,ac(iYr)])
    
    print(plot(mp))
    save(mp,jb,file=file.path(path,paste("mp",ftar,iYr,".RData",sep="_")))
    
    ## HCR
    par=hcrParam(ftar =ftar*refpts( mp)["fmsy"],
                 btrig=btrig*refpts(mp)["bmsy"],
                 fmin =fmin*refpts( mp)["fmsy"],
                 blim =blim*refpts( mp)["bmsy"])
    ref=mean(aaply(catch(om)[,ac(iYr-seq(interval))],2,sum))
    tac=hcr(mp,refs=par,hcrYrs=iYr+seq(interval),tac=TRUE)
    tac[is.na(tac)]=1
    tac[]=qmax(tac,ref*bndTac[1])
    tac[]=qmin(tac,ref*bndTac[2])

    #### OM Projectionfor TAC
    om =fwd(om,catch=tac,sr=eq,residuals=exp(sr_deviances),maxF=maxF)
    print(plot(window(om,end=iYr+3)))
    }
  
  cat("==\n")
  
  return(om)}

