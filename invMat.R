#naive matrix inversion

invMat<-function(X,eps=1e-12){
  eig.X=eigen(X,symmetric = T)
  P<- eig.X[[2]] 
  l<- eig.X[[1]] 
  ind<- l>eps
  l[ind]<- 1/l[ind] 
  l[!ind]<- 0
  ans<-P%*%diag(l,nrow=length(l))%*%t(P)
  return(ans) 
}
