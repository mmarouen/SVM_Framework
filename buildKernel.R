#compute kernel matrix

buildKernel<-function(ker, #kernel name: 'linear','gaussian','polynomial'
                      X, #input data
                      Z, # 2nd matrix
                      y, #response vector
                      degree=3, #poly degree
                      gamm=0.1, #gaussian parameter
                      pred=FALSE #binary used for prediction
                      ){
  H=c()
  S=c()
  if(ker=='linear') degree=1
  if(ker%in%c('polynomial','linear')) kermat=((Z%*%t(X)+1)^degree)-1
  if(ker=='gaussian') kermat=apply(X,1,function(x) exp(-gamm*rowSums(t(t(Z)-x)^2)))
  if(!pred){
    H=cbind(1,kermat)
    S=cbind(0,kermat)
    S=rbind(0,S)
  }
  return(list(H=H,KerMat=kermat,S=S))
}
