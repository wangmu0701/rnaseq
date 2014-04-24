




# poission distribution
p.fun=function(lambda,y){
#    print(exp(-lambda)*lambda^y)
	return( exp(-lambda)*lambda^y )
}


# Q function ( intensity part )
f.ind=function(beta,abc,start,end,exon,y){
l=matrix(0,nrow=NT,ncol=2)
l[,2]=beta[1]; l[,1]=beta[2:(1+NT)]
if(length(abc)==2){
	a=abc[1];b=abc[2]
	tmp=0
	for(i in 1:2){
		for(j in 1:2){
#            print(get(paste("Eu",i,j,exon,sep="")))
		tmp=tmp+sum( get(paste("Eu",i,j,exon,sep="")) * (-l[a,i]-l[b,j]+ y[start:end]*log(l[a,i]+l[b,j]) ) )
		}
	}
}
if(length(abc)==1){
#printf("never see this")
	a=abc[1]
	tmp=0
	for(i in 1:2){
		tmp=tmp+sum( get(paste("Eu",i,exon,sep="")) * (-l[a,i]+ y[start:end]*log(l[a,i]) ) )
	}
}
#print(tmp)
return(-tmp)
}


# Q function ( proportion part )
pp.ind=function(beta,abc,start,end,exon){
pi=matrix(0,nrow=NT,ncol=2)
pi[,1]=beta[1:NT]; pi[,2]=1-pi[,1]
if(length(abc)==2){
a=abc[1];b=abc[2]
tmp=0
	for(i in 1:2){
		for(j in 1:2){
			tmp=tmp+sum( get(paste("Eu",i,j,exon,sep="")) * (log( pi[a,i]*pi[b,j] )     )     ) 	
		}
	}
}
if(length(abc)==1){
	a=abc[1]
	tmp=0
	for(i in 1:2){
		tmp=tmp+sum( get(paste("Eu",i,exon,sep="")) * (log( pi[a,i] )     )   ) 	
	}
}
return(-tmp)
}



initialization=function(range){
l=pi=matrix(0,nrow=NT,ncol=2)
#l[,1]=runif(NT,min(range),max(range))
#l[,2]=runif(1,0.2,1)
#pi[,1]=runif(NT,0.2,1)
l[1,1]=7.0;l[1,2]=0.05;
l[2,1]=8.0;l[2,2]=0.05;
l[3,1]=9.0;l[3,2]=0.05;
pi[1,1]=0.4;
pi[2,1]=0.3;
pi[3,1]=0.2;
pi[,2]=1-pi[,1]
return(list(l=l,pi=pi))
}



# update Eu.  global assignment
update.Eu=function(l,pi,abc,start,end,exon,y){
#print(abc)
if(length(abc)==2){
	a=abc[1];b=abc[2]
	denom=0
	for(i in 1:2){
		for(j in 1:2){
			denom=denom+p.fun(  ( l[a,i]+l[b,j] )  , y[start:end] )   *  pi[a,i]*pi[b,j]
 #           print(denom)
		}
	}
	for(i in 1:2){
		for(j in 1:2){
			tmp=p.fun(  l[a,i]+l[b,j]  , y[start:end] )   *  pi[a,i]*pi[b,j]
			assign(paste("Eu",i,j,exon,sep=""),tmp/denom, envir = .GlobalEnv)
#            print(get(paste("Eu",i,j,exon,sep="")))
		}
	}
}
if(length(abc)==1){
printf("never see this!")
	a=abc[1]
	denom=0
	for(i in 1:2){
		denom=denom+p.fun(  ( l[a,i] )  , y[start:end] )   *  pi[a,i]
	}
	for(i in 1:2){
		tmp=p.fun(  l[a,i]  , y[start:end] )   *  pi[a,i]
		assign(paste("Eu",i,exon,sep=""),tmp/denom, envir = .GlobalEnv)
	}
}
}



ff=function(beta){
f=0
for(i in 1:N){
	if( length(trans[[i]]) >0){
		f=f+f.ind(beta,abc=trans[[i]],start=start[i], end=end[i],exon=paste("e",i,sep=""),y[i,])
	}
}
return(f)
}


pp=function(beta){
tmp=0
for(i in 1:N){
	if( length(trans[[i]]) >0){
		tmp=tmp+pp.ind(beta,abc=trans[[i]],start=start[i], end=end[i],exon=paste("e",i,sep=""))
	}
}
return(tmp)
}


like.ind=function(l,pi,abc,start,end,exon,y){
if(length(abc)==2){
	a=abc[1];b=abc[2]
	tmp=0
	for(i in 1:2){
		for(j in 1:2){
			tmp=tmp+sum( get(paste("Eu",i,j,exon,sep="")) * (log(pi[a,i]*pi[b,j])-l[a,i]-l[b,j]+ y[start:end]*log(l[a,i]+l[b,j]) ) ) 
        }
	}
}
if(length(abc)==1){
	a=abc[1]
	tmp=0
	for(i in 1:2){
		tmp=tmp+sum( get(paste("Eu",i,exon,sep="")) * (log(pi[a,i])-l[a,i]+ y[start:end]*log(l[a,i]) ) ) 					
	}
}
print(tmp);
return(tmp)
}


likelihood=function(l,pi){
tmp=0
for(i in 1:N){
	if( length(trans[[i]]) >0){
		tmp=tmp+like.ind(l,pi,abc=trans[[i]],start=start[i], end=end[i],exon=paste("e",i,sep=""),y[i,])
	}
}
return(tmp)
}





any.exon=function(y){
iter=1
change.of.para=1
eps=0.001
pre=1e-5
tmp=initialization(range=c(1,max(y)/2))
l=tmp$l;l
pi=tmp$pi;pi

print(l)
print(pi)

 l.new=matrix(0,nrow=NT,ncol=2);l.new
pi.new=matrix(0,nrow=NT,ncol=2);pi.new[,1]=pre;pi.new[,2]=1-pi.new[,1];pi.new


while(change.of.para>eps){ #print(iter);  
#print(change.of.para)

for(i in 1:N){ 
	if( length(trans[[i]]) >0){
		update.Eu(l,pi,abc=trans[[i]], start=start[i], end=end[i],exon=paste("e",i,sep=""),y[i,])
	}
}

print("likelihood=")
print(likelihood(l,pi))

fit=optim( c(l[1,2],l[,1]) , ff, NULL, method = "L-BFGS-B" ,lower=pre   )
l.new[,2]=fit$par[1]
l.new[,1]=fit$par[2:(1+NT)]

fitp=optim(as.vector(pi[,1]) , pp, NULL, method = "L-BFGS-B" ,lower=pre,upper=1-pre )
pi.new[,1]=fitp$par;pi.new[,2]=1-pi.new[,1]

change.of.para=sum(	abs(l-l.new)+ abs(pi-pi.new) 			 )

l=l.new; pi=pi.new

if(iter>500){
	 l=matrix(0,nrow=NT,ncol=2);
	pi=matrix(0,nrow=NT,ncol=2)
	break
}
iter=iter+1
}

return(list(l=as.data.frame(l),pi=as.data.frame(pi),likelihood=likelihood(l,pi)))
}


















################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################

