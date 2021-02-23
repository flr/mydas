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
  
  ts=model.frame(mcf(FLQuants(om,cpue   =function(x) apply(catch(x)/fbar(x),2,sum,na.rm=T), 
                                 catch  =function(x) apply(catch(x),2,sum,na.rm=T),
                                 ssb    =function(x) apply(ssb(  x),2,sum,na.rm=T),
                                 biomass=function(x) apply(stock(x),2,sum,na.rm=T))),drop=TRUE)
  ts$catch[is.na(ts$catch)][]=0.001
  
  ts=subset(ts,year<=maxyear)
  ts$cpue=mean(ts$biomass)*ts$cpue/mean(ts$cpue)
  if (is.null(devu)) return(ts)
  
  transform(merge(ts,devu,all.x=T),cpue=cpue*exp(residual))[,-6]}


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
                     ftar=1.0,btrig=0.7,fmin=0.01,blim=0.2,
                     start=range(om)["maxyear"]-15,end=range(om)["maxyear"],interval=3,
                     maxF=2,bndTac=c(0.5,1.5),msyCap=1.2,index="cpue",
                     path=""){ 
  
  #if ("FLQuant"%in%is(  u_deviances)) u_deviances  =FLQuants(u_deviances)
  #if ("FLQuant"%in%is(sel_deviances)) sel_deviances=FLQuants(sel_deviances)
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eq=dims(params(eq))$iter, rsdl=dims(sr_deviances)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in om")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Capacity limits maximum F
  maxF=quantile(c(apply(fbar(om),c(2),median)),prob=c(0.95))*maxF
  
  ## Loop round years
  cat("==")
  for (iYr in seq(start,range(om)["maxyear"]-interval,interval)){
    cat(iYr,", ",sep="")
    
    #### Observation Error Model, get catch and single CPUE, at the moment
    dat=oemFn(om,iYr)[,c("year","catch",index)]
    names(dat)[3]="index"
    
    if (dim(u_deviances)[4]>1){
      u=transmute(merge(dat[,c("year","index")],as.data.frame(u_deviances,drop=T),all.x=T),Yr=year,season=season,Index=index*exp(data))}
    else 
      u=transmute(merge(dat[,c("year","index")],as.data.frame(u_deviances,drop=T),all.x=T),Yr=year,season=1,     Index=index*exp(data))
  
    sa$cpue=cast(u,Yr~season,value="Index")[,seq(dim(u_deviances)[4]+1)] 
    sa$catch=transmute(dat, Yr=year,Total=catch)
    
    #### Management Procedure
    ## Assessment
    sa$se   =sa$cpue
    sa$se[,-1]=0.1
  
    jb=do.call("build_jabba",sa)
    jb=fit_jabba(jb,
                 init.values=TRUE,
                 init.K=sa$K.prior[1],
                 init.r=sa$r.prior[1],
                 init.q=rep(1,dim(sa$cpue)[2]-1),
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
    
    save(mp,jb,file=file.path(path,paste("mp",iYr,".RData",sep="_")))
    
    ## HCR
    par=hcrParam(ftar =ftar*refpts( mp)["fmsy"],
                 btrig=btrig*refpts(mp)["bmsy"],
                 fmin =fmin*refpts( mp)["fmsy"],
                 blim =blim*refpts( mp)["bmsy"])
    ref=mean(aaply(catch(om)[,ac(iYr)],2,sum))
    tac=hcr(mp,refs=par,hcrYrs=iYr+seq(interval),tac=TRUE)
    tac[is.na(tac)]=0.001
    tac[]=qmax(tac,ref*bndTac[1])
    tac[]=qmin(tac,min(ref*bndTac[2],refpts(mp)["msy"]*msyCap))
    mp=fwd(mp,catch=tac)
  
    #### OM Projectionfor TAC
    om =fwd(om,catch=tac,sr=eq,residuals=exp(sr_deviances),maxF=maxF)
    
    if (FALSE){
    ggplot(melt(refTimeSeries(mse,eq,mp,om),id=c("what","year")))+
      geom_hline(aes(yintercept=1),col="red")+
      geom_line(aes(year,value,col=what))+
      facet_grid(variable~.,scale="free")+
      #scale_color_manual("",values=c("blue","green"))+
      xlab("Year")+ylab("")+theme_bw()
      }
   }
  
  cat("==\n")
  
  return(om)}

