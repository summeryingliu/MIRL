#choose probability threshold by cross validation
#m fold cross validation
#B original data with missing
#n is number of obs. in the original data
#im is number of imputation
#p number of variables, returns best prob threshold
threshold<-function(x,y,q2,im=5,m=4,thr=c(0.5,0.6,0.7,0.8,0.9)){
  rand=sample(n)%%m+1
  MSE=matrix(0,length(thr),m)
  n=dim(x)[1]
  p=dim(x)[2]
  for (i in 1:m){
    X=x[rand!=i,]
    Y=y[rand!=i]
    R=mirl(X,Y,q2,im=5,E=NULL)
    PP=R[[1]]
    for (j in 1:length(thr)){
      if (max(PP[1:p+1])>thr[j]){
        W=(PP>thr[j])[2:(p+1)]
      	test=na.omit(cbind(x[rand==i,W],y[rand==i]))
        yt=test[,ncol(test)]
        xt=cbind(rep(1,nrow(test)),test[,-ncol(test)])
        cotest=R[[2]][PP>thr[j]]
        MSE[j,i]=apply((yt-xt%*%cotest)^2,2,mean)
        #NT[j,i]=d[1]
      }
      else{
        yt=y[rand==i]
        #NT[j,i]=n/m-sum(is.na(yt))
        MSE[j,i]=var(yt)*(n-1)/n
      }
    }}
  mimi=apply(MSE,1,mean,na.rm=1)
  #se=sqrt(apply(MSE,1,var)/m)
  best=which.min(mimi)
  s=best
  best=thr[s]
}
