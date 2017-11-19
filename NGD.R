#Newton-raphson optimization
NGD<-function(X, #input data
              resp, #response matrix
              Ker, #kernel matrix
              classifier, #'SVM' or 'Softmax' classifier
              lambd, #regularizer=1/C
              tol, #tolerence
              epochs #n_iterations
              ){
  N=nrow(resp)
  K=ncol(resp)
  sv=rep(TRUE,N)
   if(N>=1000){
     N1=N%/%2
     X1=X[1:N1,]
     resp1=as.matrix(resp[1:N1,])
     Ker1=Ker[1:N1,1:N1]
     out=NGD(X1,resp1,Ker1,classifier,lambd,tol,epochs)
     sv[1:N1]=out$sv
  }
  betahat=c()
  objCurve=c()
  t=1
  if(classifier=="SVM" & K==1){
    repeat{
      betahat=rep(0,N)
      intercept=0
      Ksv=Ker[sv,sv]
      Ysv=c(0,resp[sv,])
      IntMat=Ksv+lambd*diag(sum(sv))
      IntMat=rbind(1,IntMat)
      IntMat=cbind(1,IntMat)
      IntMat[1,1]=0
      bsv=solve(IntMat)%*%Ysv
      intercept=bsv[1]
      betahat[sv]=bsv[-1]
      yhat=Ker%*%betahat+intercept
      obj=objective(classifier,betahat =betahat,kerMat = Ker,y = resp,f_x = yhat,lambd = lambd)
      objCurve=c(objCurve,obj)
      sv_new=(resp*yhat) <1
      if(all(sv_new==sv)|t>50) break
      sv=sv_new
      print(paste(round(as.numeric(obj),3),sum(sv)))
      t=t+1
    }
    betahat=c(intercept,betahat)
  }
  if(classifier=="Softmax" & K==1){
    obj1=0
    repeat{
      betahat_old=rep(0,N)
      beta0_old=0
      f_x=Ker%*%betahat_old+beta0_old
      yhat=softmax(f_x)
      obj=objective(classifier,betahat =betahat_old,kerMat = Ker,y = resp,f_x = yhat,lambd = lambd)
      objCurve=c(objCurve,obj)
      grad_b=Ker%*%(lambd*betahat_old+yhat-resp)
      grad_b0=sum(yhat-resp)
      grad=c(grad_b0,grad_b)
      hessian_b=Ker%*%(lambd*diag(N)+diag(c(yhat*(1-yhat)))%*%Ker)
      hessian_b0=sum(yhat*(1-yhat))
      hessian_bb0=Ker%*%(yhat*(1-yhat))
      hessian=matrix(ncol=N+1,nrow=N+1)
      hessian[2:(N+1),2:(N+1)]=hessian_b
      hessian[1,1]=hessian_b0
      hessian[2:(N+1),1]=hessian_bb0
      hessian[1,2:(N+1)]=t(hessian_bb0)
      betahat=c(beta0_old,betahat_old)-invMat(hessian,eps = 1e-7)%*%grad
      thd=abs(obj1-obj)
      # if(t%%5==0){
      #   print(paste(t,obj))
      #   print(thd)
      # }
      if(thd<tol) break
      # if(t==6) break
      betahat_old=betahat[-1]
      beta0_old=betahat[1]
      t=t+1
      obj1=obj
    }
    sv_new=0
  }
  
  return(list(betahat=betahat,n_iter=t,gradsErr=c(),lossCurve=objCurve,sv=sv_new))
}
