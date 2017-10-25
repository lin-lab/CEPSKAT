Get_Betas <-
function(MAF,MAF.cutoff,prop.caus,c,prop.pos){
  p = length(MAF)
  Beta = rep(0,p)
  for(i in 1:p){
    if(MAF[i]<MAF.cutoff && rnorm(1)<prop.caus && MAF[i]>0){
      Beta[i] = -c*log10(MAF[i])*(2*rbinom(1,1,prop.pos)-1)
    }
  }
  return(Beta)
}
