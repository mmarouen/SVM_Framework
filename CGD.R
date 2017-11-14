#conjugate gradient descent classifier

CGD<-function(X, #input matrix
              resp, #response matrix
              Ker,H,S, #kernel matrices
              classifier, #'SVM' or 'Softmax' classifier
              lambd, #regularizer
              tol, #tolerence value
              epochs #n_iterations
              ){
  N=nrow(resp)
  K=ncol(resp)
  betahat=c()
  t=1
  loss=c()
  if(classifier=="SVM" & K==1){
    sv1=rep(TRUE,N)
    I0=diag(sv1)
    A=2*(lambd*S+t(H)%*%I0%*%H)
    b=2*t(H)%*%I0%*%resp
    beta_old=rep(0,(N+1))
    yhat=H%*%beta_old
    obj0=objective(classifier,beta_old,Ker,resp,yhat,lambd)
    loss=obj0
    repeat{
      r=b-A%*%beta_old
      alpha=as.numeric(t(r)%*%r/(t(r)%*%A%*%r))
      betahat=beta_old+alpha*r
      yhat=H%*%betahat
      sv=resp*yhat<1
      I0=I0[sv,sv]
      A=2*(lambd*S+t(H)%*%I0%*%H)
      b=2*t(H)%*%I0%*%resp
      obj1=objective(classifier,betahat,Ker,resp,yhat,lambd)
      loss=c(loss,obj1)
      # if(all(sv1==sv)) break
      if(abs(obj1-obj0)) break
      t=t+1
      beta_old=betahat
      sv1=sv
      obj0=obj1
    }
  }
  if(classifier=="Softmax"){}
  return(list(betahat=betahat,n_iter=t,gradsErr=c(),lossCurve=loss,sv=sv))
}
