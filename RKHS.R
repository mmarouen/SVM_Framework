#Gateway function to implement solve svm

RKHS<-function(Input, #covariates
              Y, #response vector
              opType='Classification', #classification or regression
              classifier='SVM', #can use SVM or softmax classifiers
              kernel='linear', #kernel type
              degree=3, #polynomial degree
              gamm=1, #gamma used in gaussian kernel
              C=0.05, #cost parameter
              learning_rate=NULL, #learning rate
              tol=1e-8, #tolerence value
              epochs=400, #number of iterations
              traceObj=F, #trace loss function for gradient descent
              gradCheck=F, #analytical gradient check
              optMode='NGD' #NGD=Newton-raphson optimization, CGD:Conjugate gradient descent
              ){
    
  rsp=transformResponse(Y,opType)
  resp=rsp$respMat
  X=as.matrix(Input)
  yhat=c()
  betahat=c()
  betahatMat=matrix()
  lossCurve=c()
  KerMat=buildKernel(kernel,X,X,resp,degree=degree,gamm=gamm)
  sv=c()
  
  Ker=KerMat$KerMat
  H=KerMat$H
  S=KerMat$S
  Ker=KerMat$KerMat
  
  lambd=1/C

  if(opType=='Regression' & classifier=='LS'){
    inv=invMat(t(H)%*%H+lambd*S)
    betahat=inv%*%t(H)%*%resp
  }
  
  if(opType=='Classification'){
    if(optMode=='BGD') out=BGD(X,resp,H,S,classifier,lambd,learning_rate,tol,epochs,traceObj,gradCheck)
    if(optMode=='SGD') out=SGD(X,resp,H,S,classifier,lambd,tol,epochs,traceObj)
    if(optMode=='NGD') out=NGD(X,resp,Ker,classifier,lambd,tol,epochs)
    if(optMode=='CGD') out=CGD(X,resp,Ker,H,S,classifier,lambd,tol,epochs)
    
    betahat=out$betahat
    lossCurve=out$lossCurve
    n_iter=out$n_iter
    gradsErr=out$gradsErr
  }
  
  f_x=H%*%betahat
  sv=out$sv
  out=transformOutput(f_x,opType,rsp$CL)
  yhat=out$yhat
  modelargs=list(kername=kernel,gamm=gamm,degree=degree,X=X,y=resp,epochs=epochs,
                 lambd=lambd,classes=rsp$CL,opType=opType)
  return(list(yhat=yhat,fx=f_x,yhatMat=out$yhatMat,beta=betahat,sv=sv,rkhsargs=modelargs,
              lossCurve=lossCurve,n_iter=n_iter,gradsErr=gradsErr))
}
