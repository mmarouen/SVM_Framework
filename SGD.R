#stochastic gradient descent

SGD<-function(X,y,H,S,ll,lambd=0.05,tol=1e-8,ep=200,traceObj=F){
  K=ncol(y)
  N=nrow(X)
  lossCurve=c()
  n_iter=1
  gradsErr=c()
  if(K==1) betahat=matrix(0,nrow=N+1,ncol=1)
  if(K>1) betahat=matrix(0,ncol=K,nrow=(N+1))
  for(t in 1:ep){
    i=sample(1:N,1)
    lr_t=1/(lambd*t)
    fx_i=H[i,]%*%betahat
    grad=lambd*S%*%betahat
    if(K==1 & (y[i,1]*fx_i<1) & ll=='hingeLoss') grad=grad-y[i,1]*H[i,]
    #if(K==1 & (lr_t*y[i,1]*fx_t<1) & ll=='hingeLoss') betahat[i,1]=betahat[i,1]+1
    if(K>1 & sum(fx_i[-which(y[i,]==1)]+1>fx_i[which(y[i,]==1)])>0 & ll=='hingeLoss'){
      k=which(y[i,]==1)
      s=sum(fx_i[-k]+1>fx_i[k])
      grad[,-k]=grad[,-k]+matrix(rep(H[i,],(K-1)),ncol=K-1)
      grad[,k]=grad[,k]-s*H[i,]
    }
    betahat=betahat-lr_t*grad
    
  }
  return(list(beta=betahat,lossCurve=lossCurve,n_iter=n_iter,gradsErr=gradsErr))
}
