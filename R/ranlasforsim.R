#xy is original n*p data with missing
#E is mice generated data
#col 1:p is x, p+1 is y
#p is number of variables
#im number of imputation
mirl<-function(x=NULL,y,q2,im=5,E=NULL,lam=exp(seq(from=log(0.55),to=log(0.001),length.out=70))){
if (!is.mids(E)) {E=mice(x,m=im,seed=123457,printFlag=FALSE)}
p=dim(x)[2]
n=dim(x)[1]
M=matrix(0,10*im,p+1)
y=as.matrix(y)
#for sth imputation
for (s in 1:im)
{MM=complete(E,s)
 x=as.matrix(MM)
 nn=length(y)
 one <- rep(1, nn)
 meanx <- drop(one %*% x)/nn
 xc <- scale(x, meanx, FALSE)         # first subtracts meanrep
 normx <- sqrt(drop(one %*% (xc^2)))
 names(normx) <- NULL
 xs <- scale(xc, FALSE, normx)
for (j in 1:10)
{  w=sample.int(length(y),replace =TRUE)
   x1=xs[w,]
   y1=y[w]
   fit1=cv.glmnet(x1,y1)
   cod=as.matrix(predict(fit1,s=fit1$lambda.1se,type="coeff"))
   if((sum(cod!=0)==1)|(sum(cod!=0)>n)){ M[(s-1)*10+j,1]=cod[1]}
   if(sum(cod!=0)!=1){M[(s-1)*10+j,which(cod!=0)]=coef(lm(y1~x1[,which(cod!=0)[-1]-1]))}
}
}
#M take records of coef for im imputation and 10 bootstrap samples to compute importance measure
m=dim(M)[2]
Impms=as.matrix(colSums(abs(M[,2:m]))/(10*im))
#print(Impms)
N=array(0,dim=c(50*im,p+1,length(lam)))
NN=matrix(0,50*im,p+1)

#N take records of coef for stability selection
#NN for compute coefficient
for (s in 1:im)
{
 MM=complete(E,s)
 x=as.matrix(MM[,1:p])
for (j in 1:50)
{  w=sample.int(as.integer(length(y)/2),replace =TRUE)
   v=sample.int(p,size=q2,prob=Impms)
   x1=x[w,v]
   y1=y[w]
   fit1=cv.glmnet(x1,y1,lambda=lam)
    N[(s-1)*50+j,c(1,v+1),1:dim(predict(fit1,s=fit1$lambda,type="coeff"))[2]]=as.matrix(predict(fit1,s=fit1$lambda,type="coeff"))

   #N is used to compute stability selection probability
   bob=as.matrix(predict(fit1,s=fit1$lambda.1se,type="coeff"))
   if((sum(bob!=0)==1)|(sum(bob!=0)>n)){ NN[(s-1)*50+j,1]=bob[1]}
   else{NN[(s-1)*50+j,c(1,v+1)][which(bob!=0)]=coef(lm(y1~x1[,which(bob!=0)[-1]-1]))
   }#cbind(NN[(s-1)*50+j,c(1,v+1)],bob)
}
}
P=apply(abs(sign(N)),c(2,3),sum)/(50*im)
Probability=apply(P,1,max)
coef=colSums(NN)/(50*im)
r=list(Probability=Probability,coef=coef)
}
#PP is stablity selection probability
#co is the coefficient estimated

