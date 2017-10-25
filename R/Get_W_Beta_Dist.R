Get_W_Beta_Dist <-
function(MAF,weights.beta){
  return(diag(dbeta(MAF,weights.beta[1],weights.beta[2])^2))
}
