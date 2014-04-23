# simulation



generate.mpi=function(p1,k1,lam1,lam2){
	tmp=rmultinom(k1,size=1,prob=c(p1,1-p1))
	intensity1=rep(lam2,k1)
      intensity1[as.numeric(tmp[1,])==1]=lam1
	return(  rpois(k1,lambda=intensity1) )
}


# generate reads
generate.filter.reads=function(start,end,l1,l2,p1){
y.tmp=generate.mpi(p1,n,l1,l2)
y<-numeric(n)
for(i in start:end){
	if(y.tmp[i]==0)
		y[i]=0
	else
		y[i]=rbinom(1,size=y.tmp[i],prob=1)
}
return(y)
}



generate.data=function(intensity, proportion){
# simulate data
l=pi=matrix(0,nrow=3,ncol=2)
l[,2]=intensity[1];  l[,1]=intensity[2:4]
pi[,1]=proportion ;  pi[,2]=1-pi[,1]
y=matrix(0,nrow=2,ncol=n)
y1 <- generate.filter.reads(1,k1,     l[1,1],l[1,2],pi[1,1])
y2 <- generate.filter.reads((k1+1), n,l[2,1],l[2,2],pi[2,1])
y3 <- generate.filter.reads(1,n,      l[3,1],l[3,2],pi[3,1])
y[1,1:k1]= y1[1:k1]+y3[1:k1]
y[2,(k1+1):n]= y2[(k1+1):n]+y3[(k1+1):n]
return(y=y)
}




res=list()

for(ii in 1:100){
res[[ii]]=list()
l0=0.04
l1=8
l2=10
l3=12
p1=0.3
p2=0.2
p3=0.1
k1=2000; n=4000


#simulate 1  non-junction 
set.seed(ii)

y=generate.data(intensity=c(l0,l1,l2,l3),proportion=c(p1,p2,p3))
write.table(y, file='y.dat', row.names=FALSE, col.names=FALSE) 
write.matrix(y,file = "y.dat", sep = " ", blocksize)
stop()
apply(y,1,sum)

par(mfrow=c(1,2),mar=c(2,2,2,2))
for(i in 1:2){
	plot(y[i,],pch=19,cex=0.9,ylim=range(y))
	abline(v=k1,col="red");	abline(v=n,col="red")
}


N=dim(y)[1]; solve=1:3; NT=length(solve); solve; NT; N
type=c("1-1","2-2");type	

trans=list()
trans[[1]]=c(1,3); trans[[2]]=c(2,3); trans
start=c(1,k1+1);  end=c(k1,n) ;start;end


	# multiple intial values
	N.iter=10
	intensity=list(0)
	mixing=list(0)
	loglikelihood=numeric(N.iter)
	success=1
	while(success <= N.iter){print(paste(ii,success,sep="-"))
		tmp=any.exon(y)
		intensity[[success]]=tmp$l
		mixing[[success]]=tmp$p
		loglikelihood[success]=tmp$likelihood
		success=success+1
	}

	ind=which.max( loglikelihood );ind
	print("estimated intensities are")
    print(intensity[[ind]])
	res[[ii]][[1]]= intensity[[ind]]
	print("estimated mixing are")
    print(mixing[[ind]])
	res[[ii]][[2]]= mixing[[ind]]


}


#save(res,file="result_two_exon_three_transcripts.Rdata")
#load("result_two_exon_three_transcripts.Rdata")


lam1=numeric(0)
lam2=numeric(0)
lam3=numeric(0)

for(i in 1:19){ 
	lam1=c(lam1,res[[i]][[1]][1,1])
	lam2=c(lam2,res[[i]][[1]][2,1])
	lam3=c(lam3,res[[i]][[1]][3,1])
}

mean(lam1);mean(lam2);mean(lam3)
sd(lam1);sd(lam2);sd(lam3)


par(mfrow=c(1,3),mar=c(2.9,4,1,1))
plot(lam1,cex=0.8,pch=19,col="red",ylab="lambda_12",ylim=c(4,19))
plot(lam2,cex=0.8,pch=19,col="green",ylab="lambda_22",ylim=c(4,19))
plot(lam3,cex=0.8,pch=19,col="blue",ylab="lambda_32",ylim=c(4,19))




