#performs SGD using pegasos implementation

pegasos<-function(Input, #covariates
                  y, #response vector
                  lambd, #regularization
                  epochs #number of iterations
                  ){
  p=ncol(Input)
  X=as.matrix(cbind(1,Input))
  N=nrow(X)
  betahat=rep(0,p+1)
  yhat=c()
  for(t in 1:epochs){
    i=sample(1:N,1)
    lr=1/(lambd*t)
    fx_i=sum(X[i,]*betahat)
    grad=lambd*c(0,betahat[-1])
    if(y[i]*fx_i<1) grad=grad-y[i]*X[i,]
    betahat=betahat-lr*grad
  }
  fx=X%*%betahat
  yhat=rep(-1,N)
  yhat[fx>0]=1
  return(list(yhat=yhat,betahat=betahat))
}
