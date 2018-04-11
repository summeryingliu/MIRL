\name{mirl}
\alias{mirl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mirl
}
\description{
This function produce the stability selection probability and the estimated coefficients
}
\usage{
mirl(x=NULL,y,q2,im=5,E=NULL,lam=exp(seq(from=log(0.55),to=log(0.001),length.out=70)))
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
  \item{E}{
 You can use the 'mice' function in the mice package to generate this E which is a 'mids' data type, if E is entered, E will be used and x will be ignored.
}
  \item{lam}{
  The vector of tunning parameter for each lasso implemented.
  }
}

\value{
\item{Probabilty}{This is the selection probability for each covariate. The larger the probability, the more significant the variable is related to the outcome.
Notice that the probability and coef are p+1 vectors and the first coef is the intercept term, where the probability is always zero.}
\item{coef }{the coefficient estimated}
}
\references{
Liu Y, Wang Y, Feng Y, Wall MM. VARIABLE SELECTION AND PREDICTION WITH INCOMPLETE HIGH-DIMENSIONAL DATA. The annals of applied statistics. 2016;10(1):418-450.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4872715/
}
\author{
Ying Liu
}

\examples{
#This example is a similar simulation setting in the reference paper.
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

\donttest{fit<-mirl(B,Y,p/2,im=5)
cbind(fit$coef,fit$Probability)}
}
