allplot<-
function(kk3, a=1, b=NULL,limit=NULL,xl="xl",yl="yl",lw=1.5,m=NULL,...){
kl<-dim(kk3)[2]
if(is.null(b)) b=dim(kk3)[1]
if(is.null(limit)){
plot(kk3[a:b,1],type="l",lty=1,ylab=yl,xlab=xl,lwd=lw,main=m,ylim=c(min(kk3[a:b,]),max(kk3[a:b,])),...)
o=1
for(i in 2:kl){
o=o+1
lines(kk3[a:b,i],lty=o,col=o)
}}
else{
plot(kk3[a:b,1],type="l",lty=1,ylab=yl,xlab=xl,lwd=lw,main=m,ylim=limit,...)
o=1
for(i in 2:kl){
o=o+1
lines(kk3[a:b,i],lty=o,col=o)
}
}}

#############################

aplot<-function(D,b=NULL,lamda,w=c(1,2),lt=1,co=1,xl="xl",yl="yl",m=NULL,...){
limit=c(min(D),max(D))
D<-as.matrix(D);N=dim(D)[1]
a1<-colMeans(abs(D))
a1<-a1/max(a1)
lw<-rep(w[1],N);lw[b]=w[2]
plot(D[1,]~log(lamda),type="l",col=co,lty=1,lwd=lw[1],ylab=yl,xlab=xl,main=m,ylim=limit,...) 
for(i in 2:N){
lt=lt+1;co=co+1
lines(D[i,]~log(lamda),col=co,lty=1,lwd=lw[i])}
}
####
aplot1<-function(D,b=c(3,4,5,6,10),w=c(1,2),lt=1,co=1,xl="xl",yl="yl",m=NULL,...){
limit=c(min(D),max(D))
D<-as.matrix(D);N=dim(D)[1]
a1<-colMeans(abs(D))
a1<-a1/max(a1)
lw<-rep(w[1],N);lw[b]=w[2]
plot(D[1,]~a1,type="l",col=co,lty=lt,lwd=lw[1],ylab=yl,xlab=xl,main=m,ylim=limit,...) 
for(i in 2:N){
lt=lt+1;co=co+1
lines(D[i,]~a1,col=co,lty=lt,lwd=lw[i])}
}
 
###
aplot2<-function(D,b=NULL,w=c(1,2),color="Yes",id=NULL,type="Yes",limit=NULL,xl="xl",yl="yl",m=NULL,...){
if(is.null(limit)) limit=c(min(D),max(D))
D<-as.matrix(D);N=dim(D)[1]
a1<-colMeans(abs(D))
a1<-a1/max(a1)
lw<-rep(w[1],N);  lw[b]=w[2]
co=switch(color, Yes=1:N,No=rep(1,N))
lt=switch(type,  Yes=1:N,No=rep(1,N))
plot(D[1,]~a1,type="l",col=co[1],lty=lt[1],lwd=lw[1],ylab=yl,xlab=xl,main=m,ylim=limit,...) 
abline(v=a1[id],lty=2)
for(i in 2:N){
lines(D[i,]~a1,col=co[i],lty=lt[i],lwd=lw[i])}
}


