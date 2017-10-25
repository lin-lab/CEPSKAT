SKAT_TruncatedNormalCEP_Get_Moment <-
function(c1, c2, sigma, Xa){
  a1<-(c1-Xa)/sigma
  a2<-(c2-Xa)/sigma  
  f1<-dnorm(a1)
  f2<-dnorm(a2)
  F1<-pnorm(a1)
  F2<-pnorm(a2)
  f12<-f2-f1
  F12<-F2 + 1 - F1
  M<- f12/F12 * sigma 
  V<- (a2*f2-a1*f1)/F12 + M^2 /sigma^2
  W<- (1 - V) * sigma^2
  Mu<- Xa - M
  return(list(M=M, V=V, W=W, Mu=Mu))
}
