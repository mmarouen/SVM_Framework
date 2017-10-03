#computes loss function

objective<-function(loss, #loss function computing
                    betahat, #coefficients
                    kerMat, #kernel matrix
                    y, #response vector
                    f_x, #scores matrix
                    lambd #regulrization paramter=1/C
                    ){
  N=nrow(y)
  K=ncol(y)
  obj=0
  penalty=0
  if(K==1){
    if(loss=='hingeLoss'){
      idx=1>as.vector(y)*f_x
      obj=sum(1-y[idx,]*f_x[idx,])
      penalty=t(betahat)%*%kerMat%*%betahat
    }
  }
  if(K>=3){
    if(loss=='hingeLoss'){
      for(k in 1:K){
        idx=y[,k]==1
        vec=as.vector(f_x[idx,k])
        obj=obj+sum((f_x[idx,-k]+1>vec)*(f_x[idx,-k]+1-vec))
      }
      penalty=sum(diag(t(betahat)%*%kerMat%*%betahat))
    }
  }
  obj=(1/N)*obj+(lambd/2)*penalty
  #obj=(1/N)*obj+(lambd/2)*sum(betahat[-1,]^2)
  return(obj)
}
