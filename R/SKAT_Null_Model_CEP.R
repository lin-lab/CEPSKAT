SKAT_Null_Model_CEP <-
function(formula, c1, c2, delta=0.001, gamma0="default", MAXITERNUM=20 , data=NULL, n.Resampling=0){
  
  
  out<-SKAT_TruncatedNormalCEPSolve(formula, c1, c2, delta, gamma0, MAXITERNUM,data)
  
  res = out$Y - out$Mu
  pi_1 = out$W
  X1 = out$X1
  n<-length(res)
  id_include=1:n
  
  res.out=NULL
  type.Resampling=NULL
  
  
  if(n.Resampling > 0){
    res.out<-SKAT_Get_TruncatedNormal(length(res), out$Xa, out$sigma, c1, c2, nSet=n.Resampling) - out$Mu
    
    res.out<-t(t(res.out) - colMeans(res.out))
  }
  
  re<-list(res=res, X1=X1,res.out=res.out,out_type="D", 
           n.Resampling=n.Resampling, type.Resampling=NULL, id_include=id_include, mu=out$Mu,pi_1=pi_1, coef = out$coef, Xa=out$Xa)
  
  class(re)<-"SKAT_NULL_Model"
  
  return(re)
  
}
