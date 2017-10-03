#analytical gradient checker

gradChecker<-function(loss, #loss function
                      betahat, #coefficients
                       H, #data matrix
                       S, #penalty matrix
                       y, #response matrix
                       fx, #scores matrix
                       l, #learning rate
                       eps=1e-6 #epsilon
                       ){
  I=sample(nrow(betahat),4,replace = F)
  J=sample(ncol(betahat),2,replace = T)
  options(digits = 12)
  error=rep(0,8)
  for(j in 1:2){
    for(i in 1:4){
      betaLeft=betahat
      betaRight=betahat
      betaLeft[I[i],J[j]]=betahat[I[i],J[j]]-eps
      fx_left=H%*%betaLeft
      betaRight[I[i],J[j]]=betahat[I[i],J[j]]+eps
      fx_right=H%*%betaRight
      objLeft=objective(loss,betaLeft,S,y,fx_left,l)
      objRight=objective(loss,betaRight,S,y,fx_right,l)
      numGrad=(objRight-objLeft)/(2*eps)
      analGrad=computeGrad(loss,y,fx,H,S,l,betahat)
      analGrad_ij=analGrad[I[i],J[j]]
      error[(j-1)*4+i]=(analGrad_ij-numGrad)/max(abs(analGrad_ij),abs(numGrad))
    }
  }
  return(error)
}
