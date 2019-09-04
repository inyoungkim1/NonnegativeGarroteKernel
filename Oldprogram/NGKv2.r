######
## used for sim III, p>n case, and use Strassen Inverse approach

library(MCMCpack)

###########
sc<-function(x, a=1){
x<-as.matrix(x)
n<-dim(x)[2];m=dim(x)[1]-1
for(i in 1:n){
x[,i]<-(x[,i]-mean(x[,i]))/sqrt(var(x[,i])*m) }
return(x) }
#############
#######

###################
          DIS0<-function(x){
          return(-(outer(x,x,'-'))^2)
          }
          DIS<-function(x){
          n<-dim(x)[1]
          K<-matrix(0,n,n)
          for(i in 1:dim(x)[2]){K<-K+(outer(x[,i],x[,i],'-'))^2}
          return(-K)
          }
          g0<-function(Dt,rho){
          K<-exp(tensor(Dt,rep(rho,dim(Dt)[3]),3,1));
          list(K=K,dK=K*tensor(Dt,rep(1,dim(Dt)[3]),3,1))
          }
          g1<-function(Dt,rho){
          return(exp(tensor(Dt,rho,3,1)))
          }
          g2<-function(Dt,rho,j){
          K<-exp(tensor(Dt,rho,3,1));
          list(K=K,dK=Dt[,,j]*K)
          }
          #####################
          ##########
       DIS0_p<-function(Xi){
       return(Xi%*%t(Xi))
       }
       DIS_p<-function(X){
       return(X%*%t(X))
       }
       g0_p<-function(Dt,rho){
       d=1
       K<-tensor(Dt,rep(rho,dim(Dt)[3]),3,1)+diag(0.0001,dim(Dt)[1]);
       list(K=K^d,dK=d*K^(d-1)*tensor(Dt,rep(1,dim(Dt)[3]),3,1))
       }
       g1_p<-function(Dt,rho){
       d=1
       return((tensor(Dt,rho,3,1))^d+diag(0.0001,dim(Dt)[1]))
       }
       g2_p<-function(Dt,rho,j){
       d=1
       K<-tensor(Dt,rho,3,1)+diag(0.0001,dim(Dt)[1]);
       list(K=K^d,dK=d*K^(d-1)*Dt[,,j])
       }
       ########

       #####################################
       # main function I:  initialize function:, estimate average of xi_i, alpha, and sigma^2 and tau, 
       # using REML Newton Rahpson method. 
       # retuern lambda_0=sigma^2/tau, BIC, alpha and so on
       ##########################################

ld0_Gauss<-function(Y,X,Dt,par0,it=100,ds=c(0.0001,0.001),ii,dd=c(1,1.5,5)){
#X<-as.matrix(X0);par0<-c(0.001,1,h/2);it=100;ds=c(0.0001,0.001);ii=3;dd=c(0,1.5,5)

ild=dd[2];dld=dd[3]      #ld is the parammeter of Marquardt method, ild/dld is increase/decrease factor for ld.
N=length(Y)
I<-diag(1,N);q<-dim(Dt)[3]
par<-matrix(0,it,3);colnames(par)<-c("s2","tau","h")
par[1,]<-par0

   ns1<-c("sum(rowSums(PK1*t(PK1)))","sum(rowSums(PK1*t(PdK1)))","sum(rowSums(PdK1*t(PdK1)))")
   ns2<-c("sum(rowSums(K1*Pt))^2","sum(rowSums(K1*Pt))*sum(rowSums(dK1*Pt))","sum(rowSums(dK1*Pt))^2")
   ns3<-c("-sum(diag(PK1))+PYt%*%K1%*%PY/par[1,1]","-sum(diag(PdK1))+PYt%*%dK1%*%PY/par[1,1]")
   ns<-paste("c(",paste("(N-3)*",ns1,"-",ns2,sep="",collapse=","),")",sep="")
  
   nsd<-paste("c(",paste(ns3,sep="",collapse=","),")",sep="")
 
 
 KK<-g0(Dt,par[1,3])
 K1<-KK$K;dK1<-KK$dK
 #K1<-exp(-tensor(Dt,rep(par[1,3],q),3,1));dK1<--K1*tensor(Dt,rep(1,q),3,1)  
 
  
 dV<-I+par[1,2]*K1;
 Vi=strassenInv(dV)  #Vi=solve(dV)
 Vi1<-rowSums(Vi)
 P<-Vi-Vi1%*%t(Vi1)/sum(Vi1)  #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
 Pt<-t(P);PY<-P%*%Y;PYt<-t(PY);
 par[1,1]=t(Y)%*%PY/(N-1)

 dV<-determinant(dV)$modulus[1];    #   -determinant(dV)$modulus[1] calculates log(det(dV)) in case det(dV)~0
 #lm<--0.5*(dV+log(det(t(X)%*%Vi%*%X))+(N-1)*log(par[1,1]))
 lm<--0.5*(dV+log(det(as.matrix(sum(Vi1))))+(N-1)*log(par[1,1]))
  
 PK1<-P%*%K1; dK1<-par[1,2]*dK1;PdK1<-P%*%dK1
 It<-xpnd(eval(parse(text=ns)),2)/2/(N-1)
 ld<-dd[1]*mean(diag(It))
 dldt<-0.5*eval(parse(text=nsd))

if(dd[1]==0) {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,2:3]+c<0)){c=c/2}
  par[i,2:3]=par[i-1,2:3]+c; #par[i,fix0]=fixv0 #par[i,which(par[i,]<0)]=sz;

  KK<-g0(Dt,par[i,3])
  K1<-KK$K;dK1<-KK$dK
 # K1<-exp(-tensor(Dt,rep(par[i,3],q),3,1));dK1<--K1*tensor(Dt,rep(1,q),3,1)  
  
  dV<-I+par[i,2]*K1;
  Vi=strassenInv(dV) #Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1) #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  Pt<-t(P);PY<-P%*%Y;
  par[i,1]=t(Y)%*%PY/(N-1)
  dV<-determinant(dV)$modulus[1];
 #lm=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+(N-1)*log(par[i,1]))
  lm<--0.5*(dV+log(det(as.matrix(sum(Vi1))))+(N-1)*log(par[i,1]))
  
   PYt<-t(PY);
   PK1<-P%*%K1; dK1<-par[i,2]*dK1;PdK1<-P%*%dK1 
  
   It<-xpnd(eval(parse(text=ns)),2)/2/(N-1)
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
   }
ii=ii-1;} 

} else {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It+ld*diag(2)),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,2:3]+c<0)){c=c/2}
  par[i,2:3]=par[i-1,2:3]+c; 

  KK<-g0(Dt,par[i,3])
  K1<-KK$K;dK1<-KK$dK
#  K1<-exp(-tensor(Dt,rep(par[i,3],q),3,1));dK1<--K1*tensor(Dt,rep(1,q),3,1)
  
  dV<-I+par[i,2]*K1; 
  Vi=strassenInv(dV) #Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1) #    P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  Pt<-t(P);PY<-P%*%Y;
  par[i,1]=t(Y)%*%PY/(N-1)
  dV<-determinant(dV)$modulus[1];
# lm0=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+(N-1)*log(par[i,1]))
  lm0<--0.5*(dV+log(det(as.matrix(sum(Vi1))))+(N-1)*log(par[i,1]))
  
  if(lm0>lm){ld=ld/dld
   PYt<-t(PY);
   PK1<-P%*%K1; dK1<-par[i,2]*dK1;PdK1<-P%*%dK1
   
   It<-xpnd(eval(parse(text=ns)),2)/2/(N-1)
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))
      
   lm=lm0
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
  } else
   {ld=ld*ild;if(ld==Inf|all(c==0)){i=it+1}}
 }
ii=ii-1;} }         #print(c(i,ii))}
k<-length(which(par[,1]!=0))
par[k,2]=par[k,2]*par[k,1]
V2=par[k,1]*I+par[k,2]*K1;
Vi<-strassenInv(V2);
Vi1=rowSums(Vi) #B<-solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
B<-t(Vi1)/sum(Vi1)
beta<-B%*%Y;
#beta<-B%*%Y; YXB<-Y-X%*%beta;
A<-par[k,2]*strassenInv(par[k,2]*K1/par[k,1]+I)%*%(I-X%*%B)/par[k,1]
#alpha<-par[k,2]*solve(par[k,2]*K1/par[k,1]+I)%*%YXB/par[k,1]; 
alpha<-A%*%Y
S=X%*%B+K1%*%A
Yhat<-S%*%Y
#Yhat<-X%*%beta+K1%*%alpha

#CV<-sum((Yhat-Y)/(1-diag(S)))^2
#CV<-sum((Yhat-Y)^2)+2*par[k,1]*sum(diag(S))
RSS<-sum((Yhat-Y)^2);df=sum(diag(S))
BIC<-log(RSS)+log(N)*df/N
 
list(beta=beta,alpha=alpha,par=par[1:k,],Yhat=Yhat,BIC=BIC,lh=lm,RSS=RSS,df=df)
#list(beta=beta,alpha=alpha,par=par[1:k,],lh=lm,Itt=It,gcv=t(Yhat-Y)%*%(Yhat-Y)/(mean(diag(I-A)))^2)
}
       ########################################################
       # see explanation for REML function
       #####################################
       # main function II:  estimating function: given xi_i, estimate alpha, and sigma^2 and tau, 
       # using REML Newton Rahpson method. The explaination of function can be found 
       # in the final pathway-environmental interaction project, REML method, not p-REML
       # retuern lambda_0=sigma^2/tau, BIC, alpha and so on
       ##########################################



ld10_Gauss<-function(Y,X,Dt,rho,par0,it=100,ds=c(0.0001,0.001),ii,dd=c(1,1.5,5)){
#X<-as.matrix(X0);K1=K;it=100;ds=c(0.0001,0.001);ii=3;dd=c(0,1.5,5)
ild=dd[2];dld=dd[3] ;N=length(Y)               #ld is the parammeter of Marquardt method, ild/dld is increase/decrease factor for ld.
I<-diag(1,N);q<-dim(Dt)[3];
par<-matrix(0,it,2);colnames(par)<-c("s2","tau")
par[1,]<-par0

   ns1<-c("sum(rowSums(P*t(P)))","sum(rowSums(P*t(PK1)))","sum(rowSums(PK1*t(PK1)))")
   ns3<-c("-sum(diag(P))+PYt%*%PY","-sum(diag(PK1))+PYt%*%K1%*%PY")
   ns<-paste("c(",paste(ns1,sep="",collapse=","),")",sep="")
   nsd<-paste("c(",paste(ns3,sep="",collapse=","),")",sep="")
 
 
 K1<-g1(Dt,rho)
#K1<-exp(-tensor(Dt,rho,3,1));  
 dV<-par[1,1]*I+par[1,2]*K1;
 Vi=strassenInv(dV) #Vi=solve(dV)
 Vi1<-rowSums(Vi)
 P<-Vi-Vi1%*%t(Vi1)/sum(Vi1) #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
 PY<-P%*%Y;PYt<-t(PY);

 dV<-determinant(dV)$modulus[1];    #   -determinant(dV)$modulus[1] calculates log(det(dV)) in case det(dV)~0
#lm=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+t(Y)%*%PY)
 lm=-0.5*(dV+log(det(as.matrix(sum(Vi1))))+t(Y)%*%PY)
 
 PK1<-P%*%K1; 
 It<-xpnd(eval(parse(text=ns)),2)/2
 ld<-dd[1]*mean(diag(It))
 dldt<-0.5*eval(parse(text=nsd))

if(dd[1]==0) {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,]+c<0)){c=c/2}
  par[i,]=par[i-1,]+c; #par[i,fix0]=fixv0 #par[i,which(par[i,]<0)]=sz;

  
  dV<-par[i,1]*I+par[i,2]*K1;
  Vi=strassenInv(dV) #  Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1)
  dV<-determinant(dV)$modulus[1];
 #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  PY<-P%*%Y;
  #lm=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+t(Y)%*%PY)
  lm=-0.5*(dV+log(det(as.matrix(sum(Vi1))))+t(Y)%*%PY)
 
   PYt<-t(PY);
   PK1<-P%*%K1;  
   It<-xpnd(eval(parse(text=ns)),2)/2
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
   }
ii=ii-1;} 

} else {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It+ld*diag(2)),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,]+c<0)){c=c/2}
  par[i,]=par[i-1,]+c; 

 
  dV<-par[i,1]*I+par[i,2]*K1;
  Vi=strassenInv(dV)  #   Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1)
  dV<-determinant(dV)$modulus[1];
 #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  PY<-P%*%Y;
 #lm0=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+t(Y)%*%PY)
  lm0=-0.5*(dV+log(det(as.matrix(sum(Vi1))))+t(Y)%*%PY)
  
  if(lm0>lm){ld=ld/dld
   PYt<-t(PY);
   PK1<-P%*%K1;# dK1<-par[i,2]*dK1;PdK1<-P%*%dK1
   It<-xpnd(eval(parse(text=ns)),2)/2
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))
      
   lm=lm0
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
  } else
   {ld=ld*ild;if(ld==Inf|all(c==0)){i=it+1}}
 }
ii=ii-1;} }         #print(c(i,ii))}
k<-length(which(par[,1]!=0))

#V2=par[k,1]*I+par[k,2]*K1;
#Vi<-solve(V2)
#beta<-solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi%*%Y; YXB<-Y-X%*%beta;
#alpha<-par[k,2]*solve(par[k,2]*K1/par[k,1]+I)%*%YXB/par[k,1]; 

#list(beta=beta,alpha=alpha,par=par[1:k,],lh=lm,Itt=It)

#V2=par[k,1]*I+par[k,2]*K1;
#Vi<-solve(V2);
#B<-solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
B<-t(Vi1)/sum(Vi1)
beta<-B%*%Y;
#beta<-B%*%Y; YXB<-Y-X%*%beta;
A<-par[k,2]*strassenInv(par[k,2]*K1/par[k,1]+I)%*%(I-X%*%B)/par[k,1]
#alpha<-par[k,2]*solve(par[k,2]*K1/par[k,1]+I)%*%YXB/par[k,1]; 
alpha<-A%*%Y
S=X%*%B+K1%*%A
Yhat<-S%*%Y
RSS<-sum((Yhat-Y)^2);df=sum(diag(S))
BIC<-log(RSS)+log(N)*df/N
list(beta=beta,alpha=alpha,par=par[1:k,],Yhat=Yhat,lh=lm,Itt=It,BIC=BIC,RSS=RSS,df=df)
}
#
#########################
ld0_poly<-function(Y,X,Dt,par0,it=100,ds=c(0.0001,0.001),ii,dd=c(1,1.5,5)){
#X<-as.matrix(X0);par0<-c(0.001,1,h/2);it=100;ds=c(0.0001,0.001);ii=3;dd=c(0,1.5,5)

ild=dd[2];dld=dd[3]      #ld is the parammeter of Marquardt method, ild/dld is increase/decrease factor for ld.
N=length(Y)
I<-diag(1,N);q<-dim(Dt)[3]
par<-matrix(0,it,3);colnames(par)<-c("s2","tau","h")
par[1,]<-par0

   ns1<-c("sum(rowSums(PK1*t(PK1)))","sum(rowSums(PK1*t(PdK1)))","sum(rowSums(PdK1*t(PdK1)))")
   ns2<-c("sum(rowSums(K1*Pt))^2","sum(rowSums(K1*Pt))*sum(rowSums(dK1*Pt))","sum(rowSums(dK1*Pt))^2")
   ns3<-c("-sum(diag(PK1))+PYt%*%K1%*%PY/par[1,1]","-sum(diag(PdK1))+PYt%*%dK1%*%PY/par[1,1]")
   ns<-paste("c(",paste("(N-3)*",ns1,"-",ns2,sep="",collapse=","),")",sep="")

   nsd<-paste("c(",paste(ns3,sep="",collapse=","),")",sep="")


 KK<-g0_p(Dt,par[1,3])
 K1<-KK$K;dK1<-KK$dK
 #K1<-exp(-tensor(Dt,rep(par[1,3],q),3,1));dK1<--K1*tensor(Dt,rep(1,q),3,1)


 dV<-I+par[1,2]*K1;
 Vi=strassenInv(dV)  #Vi=solve(dV)
 Vi1<-rowSums(Vi)
 P<-Vi-Vi1%*%t(Vi1)/sum(Vi1)  #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
 Pt<-t(P);PY<-P%*%Y;PYt<-t(PY);
 par[1,1]=t(Y)%*%PY/(N-1)

 dV<-determinant(dV)$modulus[1];    #   -determinant(dV)$modulus[1] calculates log(det(dV)) in case det(dV)~0
 #lm<--0.5*(dV+log(det(t(X)%*%Vi%*%X))+(N-1)*log(par[1,1]))
 lm<--0.5*(dV+log(det(as.matrix(sum(Vi1))))+(N-1)*log(par[1,1]))

 PK1<-P%*%K1; dK1<-par[1,2]*dK1;PdK1<-P%*%dK1
 It<-xpnd(eval(parse(text=ns)),2)/2/(N-1)
 ld<-dd[1]*mean(diag(It))
 dldt<-0.5*eval(parse(text=nsd))

if(dd[1]==0) {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,2:3]+c<0)){c=c/2}
  par[i,2:3]=par[i-1,2:3]+c; #par[i,fix0]=fixv0 #par[i,which(par[i,]<0)]=sz;

  KK<-g0_p(Dt,par[i,3])
  K1<-KK$K;dK1<-KK$dK
 # K1<-exp(-tensor(Dt,rep(par[i,3],q),3,1));dK1<--K1*tensor(Dt,rep(1,q),3,1)

  dV<-I+par[i,2]*K1;
  Vi=strassenInv(dV) #Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1) #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  Pt<-t(P);PY<-P%*%Y;
  par[i,1]=t(Y)%*%PY/(N-1)
  dV<-determinant(dV)$modulus[1];
 #lm=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+(N-1)*log(par[i,1]))
  lm<--0.5*(dV+log(det(as.matrix(sum(Vi1))))+(N-1)*log(par[i,1]))

   PYt<-t(PY);
   PK1<-P%*%K1; dK1<-par[i,2]*dK1;PdK1<-P%*%dK1

   It<-xpnd(eval(parse(text=ns)),2)/2/(N-1)
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
   }
ii=ii-1;}

} else {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It+ld*diag(2)),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,2:3]+c<0)){c=c/2}
  par[i,2:3]=par[i-1,2:3]+c;

  KK<-g0_p(Dt,par[i,3])
  K1<-KK$K;dK1<-KK$dK
#  K1<-exp(-tensor(Dt,rep(par[i,3],q),3,1));dK1<--K1*tensor(Dt,rep(1,q),3,1)

  dV<-I+par[i,2]*K1;
  Vi=strassenInv(dV) #Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1) #    P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  Pt<-t(P);PY<-P%*%Y;
  par[i,1]=t(Y)%*%PY/(N-1)
  dV<-determinant(dV)$modulus[1];
# lm0=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+(N-1)*log(par[i,1]))
  lm0<--0.5*(dV+log(det(as.matrix(sum(Vi1))))+(N-1)*log(par[i,1]))

  if(lm0>lm){ld=ld/dld
   PYt<-t(PY);
   PK1<-P%*%K1; dK1<-par[i,2]*dK1;PdK1<-P%*%dK1

   It<-xpnd(eval(parse(text=ns)),2)/2/(N-1)
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))

   lm=lm0
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
  } else
   {ld=ld*ild;if(ld==Inf|all(c==0)){i=it+1}}
 }
ii=ii-1;} }         #print(c(i,ii))}
k<-length(which(par[,1]!=0))
par[k,2]=par[k,2]*par[k,1]
V2=par[k,1]*I+par[k,2]*K1;
Vi<-strassenInv(V2);
Vi1=rowSums(Vi) #B<-solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
B<-t(Vi1)/sum(Vi1)
beta<-B%*%Y;
#beta<-B%*%Y; YXB<-Y-X%*%beta;
A<-par[k,2]*strassenInv(par[k,2]*K1/par[k,1]+I)%*%(I-X%*%B)/par[k,1]
#alpha<-par[k,2]*solve(par[k,2]*K1/par[k,1]+I)%*%YXB/par[k,1];
alpha<-A%*%Y
S=X%*%B+K1%*%A
Yhat<-S%*%Y
#Yhat<-X%*%beta+K1%*%alpha

#CV<-sum((Yhat-Y)/(1-diag(S)))^2
#CV<-sum((Yhat-Y)^2)+2*par[k,1]*sum(diag(S))
RSS<-sum((Yhat-Y)^2);df=sum(diag(S))
BIC<-log(RSS)+log(N)*df/N

list(beta=beta,alpha=alpha,par=par[1:k,],Yhat=Yhat,BIC=BIC,lh=lm,RSS=RSS,df=df)
#list(beta=beta,alpha=alpha,par=par[1:k,],lh=lm,Itt=It,gcv=t(Yhat-Y)%*%(Yhat-Y)/(mean(diag(I-A)))^2)
}
       ########################################################
       # see explanation for REML function
       #####################################
       # main function II:  estimating function: given xi_i, estimate alpha, and sigma^2 and tau,
       # using REML Newton Rahpson method. The explaination of function can be found
       # in the final pathway-environmental interaction project, REML method, not p-REML
       # retuern lambda_0=sigma^2/tau, BIC, alpha and so on
       ##########################################



ld10_poly<-function(Y,X,Dt,rho,par0,it=100,ds=c(0.0001,0.001),ii,dd=c(1,1.5,5)){
#X<-as.matrix(X0);K1=K;it=100;ds=c(0.0001,0.001);ii=3;dd=c(0,1.5,5)
ild=dd[2];dld=dd[3] ;N=length(Y)               #ld is the parammeter of Marquardt method, ild/dld is increase/decrease factor for ld.
I<-diag(1,N);q<-dim(Dt)[3];
par<-matrix(0,it,2);colnames(par)<-c("s2","tau")
par[1,]<-par0

   ns1<-c("sum(rowSums(P*t(P)))","sum(rowSums(P*t(PK1)))","sum(rowSums(PK1*t(PK1)))")
   ns3<-c("-sum(diag(P))+PYt%*%PY","-sum(diag(PK1))+PYt%*%K1%*%PY")
   ns<-paste("c(",paste(ns1,sep="",collapse=","),")",sep="")
   nsd<-paste("c(",paste(ns3,sep="",collapse=","),")",sep="")


 K1<-g1_p(Dt,rho)
#K1<-exp(-tensor(Dt,rho,3,1));
 dV<-par[1,1]*I+par[1,2]*K1;
 Vi=strassenInv(dV) #Vi=solve(dV)
 Vi1<-rowSums(Vi)
 P<-Vi-Vi1%*%t(Vi1)/sum(Vi1) #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
 PY<-P%*%Y;PYt<-t(PY);

 dV<-determinant(dV)$modulus[1];    #   -determinant(dV)$modulus[1] calculates log(det(dV)) in case det(dV)~0
#lm=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+t(Y)%*%PY)
 lm=-0.5*(dV+log(det(as.matrix(sum(Vi1))))+t(Y)%*%PY)

 PK1<-P%*%K1;
 It<-xpnd(eval(parse(text=ns)),2)/2
 ld<-dd[1]*mean(diag(It))
 dldt<-0.5*eval(parse(text=nsd))

if(dd[1]==0) {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,]+c<0)){c=c/2}
  par[i,]=par[i-1,]+c; #par[i,fix0]=fixv0 #par[i,which(par[i,]<0)]=sz;


  dV<-par[i,1]*I+par[i,2]*K1;
  Vi=strassenInv(dV) #  Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1)
  dV<-determinant(dV)$modulus[1];
 #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  PY<-P%*%Y;
  #lm=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+t(Y)%*%PY)
  lm=-0.5*(dV+log(det(as.matrix(sum(Vi1))))+t(Y)%*%PY)

   PYt<-t(PY);
   PK1<-P%*%K1;
   It<-xpnd(eval(parse(text=ns)),2)/2
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
   }
ii=ii-1;}

} else {i=2;
 while(ii>0){st=TRUE                  #st is the logic indicator for stop creterion
 while((i<=it)&&st){
  c<-tryCatch(solve(It+ld*diag(2)),error=function(e) diag(0,2))%*%dldt
  while(any(par[i-1,]+c<0)){c=c/2}
  par[i,]=par[i-1,]+c;


  dV<-par[i,1]*I+par[i,2]*K1;
  Vi=strassenInv(dV)  #   Vi=solve(dV)
  Vi1<-rowSums(Vi)
  P<-Vi-Vi1%*%t(Vi1)/sum(Vi1)
  dV<-determinant(dV)$modulus[1];
 #P<-Vi-Vi%*%X%*%solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
  PY<-P%*%Y;
 #lm0=-0.5*(dV+log(det(t(X)%*%Vi%*%X))+t(Y)%*%PY)
  lm0=-0.5*(dV+log(det(as.matrix(sum(Vi1))))+t(Y)%*%PY)

  if(lm0>lm){ld=ld/dld
   PYt<-t(PY);
   PK1<-P%*%K1;# dK1<-par[i,2]*dK1;PdK1<-P%*%dK1
   It<-xpnd(eval(parse(text=ns)),2)/2
   dldt<-0.5*eval(parse(text=gsub("[1,1]","[i,1]",nsd,fixed=TRUE)))

   lm=lm0
   st<-(any(abs(par[i,]-par[i-1,])>=ds[1]*(abs(par[i-1,])+ds[2])))
   i=i+1
  } else
   {ld=ld*ild;if(ld==Inf|all(c==0)){i=it+1}}
 }
ii=ii-1;} }         #print(c(i,ii))}
k<-length(which(par[,1]!=0))

#V2=par[k,1]*I+par[k,2]*K1;
#Vi<-solve(V2)
#beta<-solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi%*%Y; YXB<-Y-X%*%beta;
#alpha<-par[k,2]*solve(par[k,2]*K1/par[k,1]+I)%*%YXB/par[k,1];

#list(beta=beta,alpha=alpha,par=par[1:k,],lh=lm,Itt=It)

#V2=par[k,1]*I+par[k,2]*K1;
#Vi<-solve(V2);
#B<-solve(t(X)%*%Vi%*%X)%*%t(X)%*%Vi;
B<-t(Vi1)/sum(Vi1)
beta<-B%*%Y;
#beta<-B%*%Y; YXB<-Y-X%*%beta;
A<-par[k,2]*strassenInv(par[k,2]*K1/par[k,1]+I)%*%(I-X%*%B)/par[k,1]
#alpha<-par[k,2]*solve(par[k,2]*K1/par[k,1]+I)%*%YXB/par[k,1];
alpha<-A%*%Y
S=X%*%B+K1%*%A
Yhat<-S%*%Y
RSS<-sum((Yhat-Y)^2);df=sum(diag(S))
BIC<-log(RSS)+log(N)*df/N
list(beta=beta,alpha=alpha,par=par[1:k,],Yhat=Yhat,lh=lm,Itt=It,BIC=BIC,RSS=RSS,df=df)
}
#


       #####################################
       # main function I:  initialize function: given xi_i, estimate alpha, and sigma^2 and tau, 
       # using REML Newton Rahpson method. The explaination of function can be found 
       # in the final pathway-environmental interaction project, p-REML method, (not used)
       # retuern lambda_0=sigma^2/tau, BIC, alpha and so on
       ##########################################

## see explanation for p-REML function

              #####################################
              ######################
              #Main function III to get xi_i (here is rho_i) xi_i is uded in paper
              ####################################

rh0_Gauss<-function(Y,X0,Dt,Z,yh,beta,alpha,s2,tau,ln,it,qu=TRUE,aa,ds=c(0.0001,0.001)){
#D=Dt;beta=fit1$beta[1];alpha=fit1$alpha;Ki<-fit1$Ki;yh=fz;labmda0=fit1$par[24,1]/fit1$par[24,2];ln=10;it=100;dd=0; qu=TRUE; ds=c(0.0001,0.001)
       
rhoo<-matrix(0,ln+1,dim(Dt)[3]); 
rho<-matrix(0,it,dim(Dt)[3])
Cp<-ck<-RSS<-BIC<-ave<-df<-rep(0,ln+1);
aler<-rep(0,aa)
labmda0=s2/tau;  n<-length(Y); I<-diag(n);
Ybar<-Y-beta-labmda0*alpha/2
labmda<-c(exp(seq(log(max(t(Z)%*%(Ybar-sum(alpha)))/n),-6,length=ln)),0)

    K<-matrix(rep(1,n*n),n,n)
    Vi<-strassenInv(s2*I+tau*K)
    Vi1<-rowSums(Vi)
    B<-t(Vi1)/sum(Vi1)
    S=X0%*%B
    A<-strassenInv(K+labmda0*I)%*%(I-S)
    Yhat<-(S+K%*%A)%*%Y
     ck[1]<-1 
     RSS[1]=sum((Y-Yhat)^2)
     df[1]<-1+sum(rowSums(K*A))
     ave[1]<-mean((Yhat-yh)^2) 
     BIC[1]<-log(RSS[1])+df[1]*log(n)/n  
     Cp[1]<-RSS[1]+2*df[1]*(RSS[1]/(n-df[1]))
 
 i=2;st0=TRUE                 #st is the logic indicator for stop creterion
 while((i<=(ln+1))&&(st0||qu)){   # rint(i);print(it);print(st)
         
    rho[1,]<-rhoo[i-1,]
    i0=2; st=TRUE 
    while((i0<=it)&&st){
    rho[i0,]<-rho[i0-1,]
    for(j in sample(1:dim(Dt)[3],dim(Dt)[3],replace=FALSE))
    { 
     KK<-g2(Dt,rho[i0,],j)
     K<-KK$K; DK<-KK$dK
     DKa<-DK%*%alpha
     #rho[i0,j]<-rho[i0,j]+(t(Ybar-K%*%alpha)%*%DKa-n*labmda[i])/(t(DKa)%*%DKa)
     rho[i0,j]<-rho[i0,j]+(sum((Ybar-K%*%alpha)*DKa)-n*labmda[i])/sum(DKa^2)
     rho[i0,j]<-(1+sign(rho[i0,j]))*rho[i0,j]/2
    };#print(rho[i,])
   st<-(any(abs(rho[i0,]-rho[i0-1,])>=ds[1]*(abs(rho[i0-1,])+ds[2])))
   i0=i0+1;
   }
   ck[i]<-i0-1
   rhoo[i,]=rho[ck[i],]
   
       Vi<-strassenInv(s2*I+tau*K)
       Vi1<-rowSums(Vi)
       B<-t(Vi1)/sum(Vi1)
       S=X0%*%B
       A<-strassenInv(K+labmda0*I)%*%(I-S)
       Yhat<-(S+K%*%A)%*%Y
         RSS[i]=sum((Y-Yhat)^2)
         df[i]<-1+sum(rowSums(K*A))
         ave[i]<-mean((Yhat-yh)^2) 
         BIC[i]<-log(RSS[i])+df[i]*log(n)/n  
         Cp[i]<-RSS[i]+2*df[i]*(RSS[i]/(n-df[i]))
   #  if((BIC[i]-BIC[i-1])/(labmda[i-1]-labmda[i])>dd) st0=FALSE
   #aler[i%%aa+1]=(BIC[i]-BIC[i-1])/(log(labmda[i-1])-log(labmda[i]))>=0.000
   aler[i%%aa+1]=any(abs(BIC[i]-BIC[i-1])<=0.01*(abs(BIC[i-1])+1))||((BIC[i]-BIC[i-1])>=0.000)
   st0<-!all(aler==TRUE)
   i=i+1
 }
if(all(aler==TRUE))  k=i-1-(aa-2)  else k=i-1 
return(cbind(lambda=labmda[1:k],ck=ck[1:k],Cp=Cp[1:k],ave=ave[1:k],BIC=BIC[1:k],RSS=RSS[1:k],df=df[1:k],rho=rhoo[1:k,]))
}

####################
rh0_poly<-function(Y,X0,Dt,Z,yh,beta,alpha,s2,tau,ln,it,qu=TRUE,aa,ds=c(0.0001,0.001)){
#D=Dt;beta=fit1$beta[1];alpha=fit1$alpha;Ki<-fit1$Ki;yh=fz;labmda0=fit1$par[24,1]/fit1$par[24,2];ln=10;it=100;dd=0; qu=TRUE; ds=c(0.0001,0.001)

rhoo<-matrix(0,ln+1,dim(Dt)[3]);
rho<-matrix(0,it,dim(Dt)[3])
Cp<-ck<-RSS<-BIC<-ave<-df<-rep(0,ln+1);
aler<-rep(0,aa)
labmda0=s2/tau;  n<-length(Y); I<-diag(n);
Ybar<-Y-beta-labmda0*alpha/2
labmda<-c(exp(seq(log(max(t(Z)%*%(Ybar-sum(alpha)))/n),-6,length=ln)),0)

    K<-matrix(rep(1,n*n),n,n)
    Vi<-strassenInv(s2*I+tau*K)
    Vi1<-rowSums(Vi)
    B<-t(Vi1)/sum(Vi1)
    S=X0%*%B
    A<-strassenInv(K+labmda0*I)%*%(I-S)
    Yhat<-(S+K%*%A)%*%Y
     ck[1]<-1
     RSS[1]=sum((Y-Yhat)^2)
     df[1]<-1+sum(rowSums(K*A))
     ave[1]<-mean((Yhat-yh)^2)
     BIC[1]<-log(RSS[1])+df[1]*log(n)/n
     Cp[1]<-RSS[1]+2*df[1]*(RSS[1]/(n-df[1]))

 i=2;st0=TRUE                 #st is the logic indicator for stop creterion
 while((i<=(ln+1))&&(st0||qu)){   # rint(i);print(it);print(st)

    rho[1,]<-rhoo[i-1,]
    i0=2; st=TRUE
    while((i0<=it)&&st){
    rho[i0,]<-rho[i0-1,]
    for(j in sample(1:dim(Dt)[3],dim(Dt)[3],replace=FALSE))
    {
     KK<-g2_p(Dt,rho[i0,],j)
     K<-KK$K; DK<-KK$dK
     DKa<-DK%*%alpha
     #rho[i0,j]<-rho[i0,j]+(t(Ybar-K%*%alpha)%*%DKa-n*labmda[i])/(t(DKa)%*%DKa)
     rho[i0,j]<-rho[i0,j]+(sum((Ybar-K%*%alpha)*DKa)-n*labmda[i])/sum(DKa^2)
     rho[i0,j]<-(1+sign(rho[i0,j]))*rho[i0,j]/2
    };#print(rho[i,])
   st<-(any(abs(rho[i0,]-rho[i0-1,])>=ds[1]*(abs(rho[i0-1,])+ds[2])))
   i0=i0+1;
   }
   ck[i]<-i0-1
   rhoo[i,]=rho[ck[i],]

       Vi<-strassenInv(s2*I+tau*K)
       Vi1<-rowSums(Vi)
       B<-t(Vi1)/sum(Vi1)
       S=X0%*%B
       A<-strassenInv(K+labmda0*I)%*%(I-S)
       Yhat<-(S+K%*%A)%*%Y
         RSS[i]=sum((Y-Yhat)^2)
         df[i]<-1+sum(rowSums(K*A))
         ave[i]<-mean((Yhat-yh)^2)
         BIC[i]<-log(RSS[i])+df[i]*log(n)/n
         Cp[i]<-RSS[i]+2*df[i]*(RSS[i]/(n-df[i]))
   #  if((BIC[i]-BIC[i-1])/(labmda[i-1]-labmda[i])>dd) st0=FALSE
   #aler[i%%aa+1]=(BIC[i]-BIC[i-1])/(log(labmda[i-1])-log(labmda[i]))>=0.000
   aler[i%%aa+1]=any(abs(BIC[i]-BIC[i-1])<=0.01*(abs(BIC[i-1])+1))||((BIC[i]-BIC[i-1])>=0.000)
   st0<-!all(aler==TRUE)
   i=i+1
 }
if(all(aler==TRUE)) k=i-1-(aa-2) else k=i-1
return(cbind(lambda=labmda[1:k],ck=ck[1:k],Cp=Cp[1:k],ave=ave[1:k],BIC=BIC[1:k],RSS=RSS[1:k],df=df[1:k],rho=rhoo[1:k,]))
}
