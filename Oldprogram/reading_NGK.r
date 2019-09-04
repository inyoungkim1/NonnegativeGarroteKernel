n<-64
#############################################################
#What user needs to change: Kernel and path
name<-"Gauss"
path=paste("/FolderName/Test_",name,"_n=",n,"_p=80.csv",sep="")
#############################################################

name<-"Gauss"
path=paste("/FolderName/Test_",name,"_n=",n,"_p=80.csv",sep="")
A<-as.matrix(read.csv(path)[,-1])
 
#RSS
RSS<-mean(A[,6]/n);RSSsd<-sd(A[,6]/n)
#SE
SE<-mean(A[,4]);SEsd<-sd(A[,4])


A[which(A!=0)]<-1

count<-colSums(A)

#False Positive rate (alpha)
FP<-mean(rowSums(A[,13:87])/(rowSums(A[,13:87])+75));
FPsd<-sd(rowSums(A[,13:87])/(rowSums(A[,13:87])+75))


#False negative rate (beta)
FN<-mean((5-rowSums(A[,8:12]))/((5-rowSums(A[,8:12]))+5));
FNsd<-sd((5-rowSums(A[,8:12]))/((5-rowSums(A[,8:12]))+5))


MS<-mean(rowSums(A[,8:87]))
MSsd<-sd(rowSums(A[,8:87]))

count
cbind(FP,FPsd,FN,FNsd,MS,MSsd,RSS,RSSsd,SE,SEsd)

#############################################################
#What user needs to change: Kernel and path
name<-"Poly"
path=paste("/FolderName/Test_",name,"_n=",n,"_p=80.csv",sep="")
#############################################################

A<-as.matrix(read.csv(path)[,-1])
 
#RSS
RSS<-mean(A[,6]/n);RSSsd<-sd(A[,6]/n)
#SE
SE<-mean(A[,4]);SEsd<-sd(A[,4])


A[which(A!=0)]<-1

count1<-colSums(A)

#False Positive rate (alpha)
FP<-mean(rowSums(A[,13:87])/(rowSums(A[,13:87])+75));
FPsd<-sd(rowSums(A[,13:87])/(rowSums(A[,13:87])+75))


#False negative rate (beta)
FN<-mean((5-rowSums(A[,8:12]))/((5-rowSums(A[,8:12]))+5));
FNsd<-sd((5-rowSums(A[,8:12]))/((5-rowSums(A[,8:12]))+5))


MS<-mean(rowSums(A[,8:87]))
MSsd<-sd(rowSums(A[,8:87]))

count1
cbind(FP,FPsd,FN,FNsd,MS,MSsd,RSS,RSSsd,SE,SEsd)


#############
 par(mfrow=c(1,1),oma=c(0,1,4,0),mar=c(8,4.5,6,2))
 plot(count[-(1:7)]/400,xlab="Predictor Index", ylab="Probability",cex=1.,cex.lab=1.5,cex.axis=1.5)
 points(count1[-(1:7)]/400,pch=20)
 legend(locator(1),cex=1.5, c("NGK Gauss","NGK linear poly"),pch=c(1,20))

 


	
