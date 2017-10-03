#svm predictor

predictRkhs<-function(Input, #test data
                      rkhs #rkhs model
                      ){
  X=as.matrix(Input)
  gamm=rkhs$rkhsargs$gamm
  degree=rkhs$rkhsargs$degree
  Xtrain=rkhs$rkhsargs$X
  y=rkhs$rkhsargs$y
  epochs=rkhs$rkhsargs$epochs
  kerName=rkhs$rkhsargs$kername
  lambd=rkhs$rkhsargs$lambd
  out=buildKernel(kerName,Xtrain,X,y,degree,gamm,pred=TRUE)
  Ker=out$KerMat
  alphahat=rkhs$beta
  beta0=alphahat[1]
  fx=beta0+Ker%*%as.matrix(alphahat[-1])
  yhat=rep(-1,nrow(Input))
  yhat[fx>0]=1
  return(yhat)
}
