SKAT_Get_TruncatedNormal <-
function(n, mu, sigma, c1, c2, nSet, nUpper = n/2){
  
  
  prob1<-1- pnorm(c1,mean=mu, sd=sigma) 
  prob2<-pnorm(c2,mean=mu, sd=sigma) 
  
  prob<-prob1/(prob1 + prob2)
  
  temp.b<-rbinom(n*nSet,1,rep(prob,nSet))# Prob x > c1
  
  lower = rep(- 10^5,n*nSet) * (1-temp.b) + c1
  upper = rep(10^5,n*nSet) * temp.b + c2
  
  re<-rtnorm(n*nSet, mean=mu, sd=sigma, lower=lower, upper=upper)
  re.M<-matrix(re, ncol=nSet, byrow=FALSE)
  
  cutoff<<-cbind(lower,upper, temp.b, mu)
  rm(re)
  rm(lower)
  rm(upper)
  
  return(re.M)
  
}
