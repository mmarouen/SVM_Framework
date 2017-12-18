#Gateway function to solve svm

RKHS<-function(Input, #covariates
              Y, #response vector
              Xtest,#test matrix
              yTest,#response test
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
    
  t1=Sys.time()
  rsp=transformResponse(Y,classifier,opType)
  resp=rsp$respMat
  X=as.matrix(Input)
  XTest=NULL
  resp1=NULL
  kerT=NULL
  if(!is.null(Xtest) & !is.null(yTest)){
    XTest=as.matrix(Xtest)
    rsp1=transformResponse(yTest,classifier,opType)
    resp1=rsp1$respMat
    kers=buildKernel(kernel,X,XTest,degree=degree,gamm=gamm)
    kerT=kers$KerMat
  }
  yhat=c()
  yhatTest=c()
  betahat=c()
  betahatMat=matrix()
  lossCurve=c()
  KerMat=buildKernel(kernel,X,X,degree=degree,gamm=gamm)
  sv=c()
  
  Ker=KerMat$KerMat
  H=KerMat$H
  S=KerMat$S
  Ker=KerMat$KerMat
  
  lambd=1/(2*C)

  if(opType=='Regression' & classifier=='LS'){
    inv=invMat(t(H)%*%H+lambd*S)
    betahat=inv%*%t(H)%*%resp
  }
  
  if(opType=='Classification'){
    if(optMode=='BGD') out=BGD(X,resp,H,S,classifier,lambd,learning_rate,tol,epochs,traceObj,gradCheck)
    if(optMode=='SGD') out=SGD(X,resp,H,S,classifier,lambd,tol,epochs,traceObj)
    if(optMode=='NGD') out=NGD(X,resp,Ker,classifier,lambd,tol,epochs)
    if(optMode=='CGD') out=CGD(X,resp,XTest,resp1,Ker,kerT,H,S,classifier,lambd,tol,epochs)
    
    betahat=out$betahat
    lossCurve=out$lossCurve
    n_iter=out$n_iter
    gradsErr=out$gradsErr
  }
  
  f_x=H%*%betahat
  sv=out$sv
  out=transformOutput(f_x,opType,rsp$CL)
  yhat=out$yhat
  t2=Sys.time()
  modelargs=list(kername=kernel,gamm=gamm,degree=degree,X=X,y=resp,epochs=epochs,
                 lambd=lambd,classes=rsp$CL,opType=opType)
  return(list(yhat=yhat,fx=f_x,yhatMat=out$yhatMat,beta=betahat,sv=sv,rkhsargs=modelargs,
              lossCurve=lossCurve,n_iter=n_iter,gradsErr=gradsErr,duration=t2-t1))
}
