#Contains useful functions for resolution

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

#svm predictor
predictRkhs<-function(Input, #input matrix
                      rkhs #model
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
  opType=rkhs$rkhsargs$opType
  CL=rkhs$rkhsargs$classes
  out=transformOutput(fx,opType,CL)
  yhat=out$yhat
  return(yhat)
}

#computes softmax transformation while avoiding numerical overflow
softmax<-function(X){
  eps=1e-15
  Eps=1-eps
  M=max(X)
  product=apply(X,2,function(x) exp(-M-log(rowSums(exp(X-M-x)))))
  product=pmax(product,eps)
  product=pmin(product,Eps)
  return(product)
}

#converts reponse back into vector format identical to input response
transformOutput<-function(f_x, #scores
                          tt="Regression", #operation mode
                          classes #classes vector
                          ){
  yhat=c()
  if(tt=="Regression") yhat=f_x
  
  if(tt=="Classification"){
    K=ncol(as.matrix(f_x))
    if(K==1){
      yhat=rep(-1,length=nrow(f_x))
      yhat[f_x>0]=1
    }
    if(K>2){
      yhat=apply(f_x,1,function(x) classes[which.max(x)])
      yhat=as.factor(yhat)
    }
  }
  return(list(yhat=yhat,yhatMat=f_x))
}

#converts response vector into matrix format
transformResponse<-function(resp, #response vector
                            tt="Regression" #operation mode
                            ){
  classes=NULL
  if(tt=="Regression"){respMat=as.matrix(resp)}
  if(tt=="Classification"){
    K=length(unique(resp))
    if(K==2){respMat=as.matrix(resp)}
    if(K>2){
      respMat=matrix(0,nrow=length(resp),ncol=K)
      classes=sort(unique(resp))
      for (i in 1:length(resp)){respMat[i,which(resp[i]==classes)]=1}
    }
  }
  return(list(respMat=respMat,CL=classes,response=resp))
}

