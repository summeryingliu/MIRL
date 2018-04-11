\name{threshold}
\alias{threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
threshold
}
\description{
This function picked the best cutting threshold for the selection probability. To pick a final set of variables. The criterion is based on cross validation, to find the threshold that minimize the prediction error in validation set. Notice that this function will run slower than the MIRL function itself since it is running mirl multiple times.

Instead the user can also self define the selection probability threshold. to save the complexity of choosing by CV.
}
\usage{
threshold(x, y, q2, im,m=4,thr=c(0.5,0.6,0.7,0.8,0.9))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
This is the n by p design matrix with missing entries.
}
  \item{y}{
  This is a n by 1 vector of the outcome, the outcome should be non missing. Please delete any sample with missing outcome.
  }
  \item{q2}{
  This is the number of variables to be bootstraped, recommended size is p/2.
}
  \item{im}{
 Number of multiple imputation, increase of im will increase time cost. Default is 5.
}

  \item{m}{
Number of folds for cross validation. The default is 4. Notice that increase this number will increase time linearly.
}
  \item{thr}{
The threshold grid to pick from.}
}

\value{
\item{best}{The best threhold.}
}

\author{
Ying Liu
}

\examples{
#This example is a simulation setting in the reference paper.
cor=0.6
prob=0.02
p=10
n=200 #sample size
Sigma=matrix(cor,p,p)#correlation of predictors
diag(Sigma)=1
mu=numeric(p)
set.seed(3)
#C is the complete design matrix without missing
C=mvrnorm(n,mu,Sigma)
#The missing indicator matrix
A<-matrix(rbinom(n*p,size=1,prob),n,p)
A[,c(1,4,6)]=0 #columns without missing
p1=inv.logit(C[,1]+C[,6]-2)
A[,5]=rbinom(n,size=1,p1) #Missing at Random
p2=inv.logit(-C[,1]-0.5*C[,6]-2)
A[,10]=rbinom(n,size=1,p2)
p3=inv.logit(C[,4]-2)
A[,9]=rbinom(n,size=1,p3)
beta=numeric(p)
beta[1:6]=c(0.1,0.2,0.5,-0.3,-.4,-0.5)*5
ct=c(0,beta)
#generating Y
Y=C\%*\%beta+rnorm(n)
B=C
B[A==1]=NA


\donttest{
best<-threshold(B,Y,q2=p/2,m=3,thr=c(0.75,0.85))
fit<-mirl(B,Y,3,p/2,im=5)
#the column number for selected variables
select=which(fit$Probability>best)
}
}
