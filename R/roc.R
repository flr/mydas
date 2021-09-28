
roc<-function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), 
             FPR=cumsum(!labels)/sum(!labels),
             labels,
             reference=sort(scores))}
# 
# 
# library(mydas)
# 
# om=rlnorm(1000)
# 
# # unbiased estimator
# mp1=om*rlnorm(1000,0.3)
# dat=mydas:::roc(om>1,mp1)
# 
# # poor estimator
# mp2=om*rlnorm(1000,0,2)
# dat2=mydas:::roc(om>1,mp2)
# 
# # biased estimator
# mp3=om*rlnorm(1000,1)
# dat3=mydas:::roc(om>1,mp3)
# 
# 
# ggplot()+
#   geom_line(aes(FPR,TPR),data=dat,col="red")+
#   geom_point(aes(FPR,TPR),data=subset(dat,(reference-1)^2==min((reference-1)^2)),size=3,,col="red")+
#   geom_line(aes(FPR,TPR),data=dat2,col="blue")+
#   geom_point(aes(FPR,TPR),data=subset(dat2,(reference-1)^2==min((reference-1)^2)),size=3,,col="blue")+
#   geom_line(aes(FPR,TPR),data=dat3,col="green")+
#   geom_point(aes(FPR,TPR),data=subset(dat3,(reference-1)^2==min((reference-1)^2)),size=3,,col="green")

# library(DescTools)
# 
# AUC(dat$FPR, dat$TPR)
# AUC(dat2$FPR,dat2$TPR)
# AUC(dat3$FPR,dat2$TPR)