#batch gradient descent using netwon optimization implementation

BGD<-function(X, #covariates
              y, #response matrix
              H, #input data
              S, #penalty data
              ll, #loss function 
              lambd=0.05, #regularization parameter lambda=1/C
              lr=NULL, #learning rate
              tol=1e-8, #tolerence
              ep=200, #number of iterations 
              traceObj=F, #trace objective function  
              gradCheck=F #gradient check
              ){
  K=ncol(y)
  N=nrow(X)
  lossCurve=c()
  n_iter=1
  gradsErr=c()
  if(K==1){ #binary classification
    betahat=matrix(0,nrow=N+1,ncol=1)
    sv=sample(c(TRUE,FALSE),N,replace = T)
    #betahat=matrix(runif(N+1),ncol=1)
    repeat{
      if(is.null(lr)) lr=1/n_iter
      
      out=updateParams(ll,y,H,S,lambd,betahat)
      betahat_new=out$betahat_new
      sv_new=out$sv_new
      # if(sum(sv_new)>1) grad=-as.matrix(colSums(as.matrix(y[sv_new,]*H[sv_new,])))
      # if(sum(sv_new)==1) grad=-y[sv_new]*H[sv_new,]
      # gradMat=(1/N)*grad+lambd*S%*%betahat
      # betahat_new=betahat-lr*gradMat
      
      if(traceObj){
        lossCurve=c(lossCurve,objective(ll,betahat,S,y,f_x,lambd))
      }
      
      if(all(sv==sv_new)) break
      # if(sqrt(sum(grad[-1,]^2))<1e-3*N*n_iter) break
      
      if(n_iter%in%c(200,800,1000) & gradCheck){
        gradsErr=rbind(gradsErr,t(gradChecker(ll,betahat,H,S,y,f_x,lambd)))
      }
      betahat=betahat_new
      sv=sv_new
      n_iter=n_iter+1
    }
  }
  if(K>=3){
    betahatMat=matrix(0,ncol=K,nrow=(N+1))
    repeat{
      if(is.null(lr)) lr=1/n_iter
      f_x=as.matrix(H%*%betahatMat)
      gradMat=computeGrad(ll,y,f_x,H,S,lambd,betahatMat)
      betahatMat_new=betahatMat-lr*gradMat
      if(traceObj){
        lossCurve=c(lossCurve,objective(ll,betahatMat,S,y,f_x,lambd))
      }
      if(n_iter>ep) break
      
      if(n_iter%in%c(200,800,1000) & gradCheck){
        gradsErr=rbind(gradsErr,t(gradChecker(ll,betahatMat,H,S,y,f_x,lambd)))
      }
      betahatMat=betahatMat_new
      n_iter=n_iter+1
    }
    betahat=betahatMat
  }
  return(list(beta=betahat,sv=sv,lossCurve=lossCurve,n_iter=n_iter,gradsErr=gradsErr))
}
