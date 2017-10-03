#Gateway function to implement solve svm

RKHS<-function(Input, #covariates
              Y, #response vector
              opType='Classification', #classification or regression
              loss='hingeLoss', #loss function used. If binary classification, then algorithm implicitely uses huberized hingeloss
              kernel='linear', #kernel type
              degree=3, #polynomial degree
              gamm=1, #gamma used in gaussian kernel
              C=0.05, #cost parameter
              learning_rate=NULL, #learning rate
              tol=1e-8, #tolerence value
              epochs=400, #number of iterations
              traceObj=F, #trace loss function for gradient descent
              gradCheck=F, #analytical gradient check
              optMode='BGD' #reslve using BGD or SGD
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

  if(opType=='Regression'){
    inv=invMat(t(H)%*%H+lambd*S)
    betahat=inv%*%t(H)%*%resp
  }
  
  if(opType=='Classification'){
    if(optMode=='BGD'){
      out=BGD(X,resp,H,S,loss,lambd,learning_rate,tol,epochs,traceObj,gradCheck)
    }
    if(optMode=='SGD'){
      out=SGD(X,resp,H,S,loss,lambd,tol,epochs,traceObj)
    }
    betahat=out$beta
    lossCurve=out$lossCurve
    n_iter=out$n_iter
    gradsErr=out$gradsErr
  }
  
  f_x=H%*%betahat
  sv=out$sv
  out=transformOutput(f_x,opType,rsp$CL)
  yhat=out$yhat
  modelargs=list(kername=kernel,gamm=gamm,degree=degree,X=X,y=resp,epochs=epochs,
                 lambd=lambd)
  return(list(yhat=yhat,fx=f_x,yhatMat=out$yhatMat,beta=betahat,sv=sv,rkhsargs=modelargs,
              lossCurve=lossCurve,n_iter=n_iter,gradsErr=gradsErr))
}
