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
