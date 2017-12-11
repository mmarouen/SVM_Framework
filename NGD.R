#Newton-raphson optimization
NGD<-function(X,resp,Ker,classifier,lambd,tol,epochs){
  N=nrow(resp)
  K=ncol(resp)
  sv=rep(TRUE,N)
  betahat=rep(0,N+1)
  if(N>=1000){
    N1=N%/%2
    X1=X[1:N1,]
    resp1=as.matrix(resp[1:N1,])
    Ker1=Ker[1:N1,1:N1]
    out=NGD(X1,resp1,Ker1,classifier,lambd,tol,epochs)
    if(classifier=="SVM") sv[1:N1]=out$sv
    if(classifier=="Softmax") betahat[1:N1]=out$betahat
  }
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

      # bsv=invMat(IntMat,eps = 1e-12)%*%Ysv
      bsv=solve(IntMat)%*%Ysv
      intercept=bsv[1]
      betahat[sv]=bsv[-1]
      yhat=Ker%*%betahat+intercept
      obj=objective(classifier,betahat =betahat,kerMat = Ker,y = resp,f_x = yhat,lambd = lambd)
      objCurve=c(objCurve,obj)
      sv_new=(resp*yhat) <1
      if(all(sv_new==sv)|t>50) break
      sv=sv_new
      # print(paste(round(as.numeric(obj),3),sum(sv)))
      t=t+1
    }
    betahat=c(intercept,betahat)
  }
  if(classifier=="Softmax" & K==1){
    obj1=0
    H=cbind(1,Ker)
    K1=cbind(0,Ker)
    K1=rbind(0,K1)
    K1[1,1]=1
    I1=cbind(0,diag(N))
    I1=rbind(0,I1)
    I2=rbind(1,diag(N))
    beta_old=betahat
    f_x=H%*%beta_old
    yhat=softmax(f_x)
    grad=2*lambd*I1%*%beta_old+I2%*%(yhat-resp)
    # thd=sum((K1%*%grad)^2)*tol
    repeat{
      grad=2*lambd*I1%*%beta_old+I2%*%(yhat-resp)
      hess=2*lambd*I1+I2%*%diag(c(yhat*(1-yhat)))%*%H
      delt=solve(hess)%*%grad
      beta_new=beta_old-delt
      f_x=H%*%beta_new
      yhat=softmax(f_x)
      t=t+1
      eps=sum(grad*delt)
      beta_old=beta_new
      if(t>15| eps<tol) break
    }
    sv_new=0
  }
  return(list(betahat=beta_new,n_iter=t,gradsErr=c(),lossCurve=objCurve,sv=sv_new))
}

