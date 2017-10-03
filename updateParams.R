#updates coefficients during a gradient descent step

updateParams<-function(loss,#loss function used
                       y, #response vector
                       H, #input matrix
                       S, #penalty matrix
                       l, #learning rate
                       betaCoef, #coefficients matrix
                       alg='newton' #optimization method
                       ){
  fx=H%*%betaCoef
  N=nrow(fx)
  K=ncol(fx)
  gradMat=matrix()
  if(K==1){
    gradMat=matrix(0,nrow=N+1,ncol=1)
    if(loss=='hingeLoss'){
      sv=1> y*fx
      I0=diag(as.vector(sv))
      betahat_new=matrix(0,nrow=N+1,ncol=1)
      M=l*diag(N)+I0%*%H[,-1]
      M=rbind(1,M)
      M=cbind(c(0,1*sv),M)
      betahat_new=solve(M)%*%c(0,sv*y)
      ll=list(betahat_new=betahat_new,sv_new=sv)
    }
    if(loss=='crossEntropy'){
      
    }
  }
  
  if(K>=3){
    gradMat=matrix(0,nrow=(N+1),ncol=K)
    if(loss=='hingeLoss'){
      for(k in 1:K){
        idx=y[,k]==1
        vec=as.vector(fx[idx,k])
        gradMat[,k]=-colSums(as.matrix(rowSums(fx[idx,-k]+1>vec)*H[idx,]))
        yk=y[!idx,-k]
        vec=as.vector(fx[!idx,k])
        gradMat[,k]=gradMat[,k]+colSums(as.matrix(rowSums((vec+1>fx[!idx,-k])*yk)*H[!idx,]))
        ll=gradMat
      }
    }
  }
  
  return(ll)
  }
