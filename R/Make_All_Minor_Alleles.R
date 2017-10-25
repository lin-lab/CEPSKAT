Make_All_Minor_Alleles <-
function(h,genotypes=FALSE,MAF){
  if(genotypes){
    IDB = which(MAF>0.5)
    h[,IDB] = 2 - h[,IDB]
    return(as.matrix(h))
  }
  else{
    IDB = which(MAF>0.5)
    h[,IDB] = 1 - h[,IDB]
    return(as.matrix(h))
  }
}
