trend<-function(x,n=12){
  if (length(x)<n) return(NA)
  
  y=x[length(x)-(n:1)+1]
  
  if (any(is.na(x))) return(NA)
  
  lm(y~x,data=data.frame(x=seq(n),y=y))$coef[2]}

twothree<-function(x,n1=2,n2=3){
  
  if (length(x)<n1+n2) return(NULL)
  
  mean(x[length(x)-seq(n1)+1])/mean(x[length(x)-(n1+n2):n2+1])}
