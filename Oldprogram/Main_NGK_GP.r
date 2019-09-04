### Gauss K and Poly K###

#####what user need to change here ####

rn<-400;N=2^6;  
setwd("/FolderName/")

outpath=getwd()

########Call subfuntions before NGK

source("R_funs.r")
source("NGKv2.r")
source("plots.r")
source("stras.r")


 
######################
i0=0
OUT<-NULL
OUT_p<-NULL

p<-80;GEN=5;       #N: sample size, p: var size, GEN:true var size; User can change
t1=0;t2=0          ## par for correlated X1 X0
X0<-rep(1,N);X0<-as.matrix(X0) ## intercept vector 1's
X<-matrix(0,N,p)
 

###############Gaussian K

 tm<-proc.time()
 while(i0<rn){
   ################
   U<-runif(N,0,1)
   for(i in 1:GEN){    X[,i]<-(runif(N,0,1)+t1*U)/(1+t1)}
   for(i in (GEN+1):p){X[,i]<-(runif(N,0,1)+t2*U)/(1+t2)}

   yh=10*cos(X[,1])+3*X[,2]^2+5*sin(X[,3])+6*exp(X[,4]/3)*X[,4]+8*cos(X[,5])+X[,5]*X[,2]*X[,1]
   
   Y<-yh+rnorm(N,0,1)
   X<-sc(X);

  
   #################
 Dt<-array(data=NA,c(dim(X)[1],dim(X)[1],dim(X)[2]))
 for(i in 1:dim(X)[2]){ Dt[,,i]<-DIS0(X[,i])}
 
 t<-DIS(X);
 h<-abs(1/(sum(t)/sum(t!=0)))
 
 fit1<-ld0_Gauss(Y,X0,Dt,par0=c(0.001,1,h/2),it=200,ds=c(0.0001,0.001),3,c(0.01,1.5,5))
 kk<-dim(fit1$par)[1]
 Z<-apply(Dt,3,function(x) x%*%fit1$alpha)
 rho<-rh0_Gauss(Y,X0,Dt,Z,yh,fit1$beta[1],fit1$alpha,fit1$par[kk,1],fit1$par[kk,2],ln=20,it=1000,qu=FALSE,aa=5,ds=c(0.0001,0.001))
 
  k<-dim(rho)[1]
  id<-which(rho[,"BIC"]==min(rho[,"BIC"]))
  OUT<-rbind(OUT,rho[id,])
  path<-paste(outpath,"/Test_Gauss_n=",N,"_p=",p,".csv",sep="")

  write.csv(OUT,path)
}#else {windows();print(er)}

proc.time()-tm




##############################Poly K

  tm<-proc.time()
 while(i0<rn){
   ################
   U<-runif(N,0,1)
   for(i in 1:GEN){    X[,i]<-(runif(N,0,1)+t1*U)/(1+t1)}
   for(i in (GEN+1):p){X[,i]<-(runif(N,0,1)+t2*U)/(1+t2)}

   yh=10*cos(X[,1])+3*X[,2]^2+5*sin(X[,3])+6*exp(X[,4]/3)*X[,4]+8*cos(X[,5])+X[,5]*X[,2]*X[,1]
   
   Y<-yh+rnorm(N,0,1)
   X<-sc(X);


 Dt<-array(data=NA,c(dim(X)[1],dim(X)[1],dim(X)[2]))
 for(i in 1:dim(X)[2]){ Dt[,,i]<-DIS0_p(X[,i])}
 
 h=2
 fit1<-ld0_poly(Y,X0,Dt,par0=c(0.001,1,h/2),it=200,ds=c(0.0001,0.001),3,c(0.01,1.5,5))
 kk<-dim(fit1$par)[1]
 Z<-apply(Dt,3,function(x) x%*%fit1$alpha)
 rho<-rh0_poly(Y,X0,Dt,Z,yh,fit1$beta[1],fit1$alpha,fit1$par[kk,1],fit1$par[kk,2],ln=20,it=1000,qu=FALSE,aa=5,ds=c(0.0001,0.001))
 
  k<-dim(rho)[1]
  id<-which(rho[,"BIC"]==min(rho[,"BIC"]))
  OUT_p<-rbind(OUT_p,rho[id,])
  i0<-dim(OUT_p)[1]
  path<-paste(outpath,"/Test_poly_n=",N,"_p=",p,".csv",sep="")

  write.csv(OUT_p,path)


  } #else {windows();print(er)}

proc.time()-tm


