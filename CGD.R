#conjugate gradient descent classifier

CGD<-function(X,resp,XT,respT, #respectively train matrix, response, test matrix, response
              Ker,kerT, #train kernel & test kernel
              H,S, #adaptations of the kernel matrix
              classifier, #classifier to be used: svm or softmax
              lambd, #regularization parameter
              tol, #tolerence factor
              epochs #max iteration
             ){
  N=nrow(resp)
  K=ncol(resp)
  betahat=c()
  t=1
  loss=c()
  
  if(classifier=="SVM" & K==1){
    sv=rep(TRUE,N)
    sv1=c()
    I0=diag(sv)
    I1=cbind(0,diag(N))
    I1=rbind(0,I1)
    I2=rbind(1,diag(N))
    A=lambd*I1+I2%*%I0%*%H
    b=I2%*%I0%*%resp
    betahat=rep(0,(N+1))
    K1=cbind(0,Ker)
    K1=rbind(0,K1)
    K1[1,1]=1
    yhat=H%*%betahat
    loss=c()
    d=2*(-A%*%betahat+b)
    g_old=d
    gnorm0=sum((K1%*%g_old)^2)
    thd=tol*gnorm0
    repeat{
      alpha1=t(d)%*%K1
      alpha=as.numeric(alpha1%*%g_old/(alpha1%*%A%*%d))
      betahat=betahat+alpha*d
      yhat=H%*%betahat
      sv1=resp*yhat<1
      I0=diag(c(sv1))
      A=lambd*I1+I2%*%I0%*%H
      b=I2%*%I0%*%resp
      g_new=-A%*%betahat+b
      numerator=as.numeric(t(g_new)%*%K1%*%g_new)
      denumerator=as.numeric(t(g_old)%*%K1%*%g_old)
      d=2*g_new+numerator*d/denumerator
      gnorm1=sum((K1%*%g_new)^2)
      if(t>50 | all(sv==sv1) | gnorm1<thd) break
      t=t+1
      sv=sv1
      g_old=g_new
    }
  }
  if(classifier=="Softmax"){}
  return(list(betahat=betahat,n_iter=t,gradsErr=c(),lossCurve=loss,sv=sv))
}
