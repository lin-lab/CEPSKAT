Sample_Haplotypes <-
function(SNP.location,Region.length,n){
  MAX.loc = max(SNP.location)
  MIN.loc = min(SNP.location)
  Start.loc = floor(runif(1)*(MAX.loc-MIN.loc-Region.length)+MIN.loc)
  End.loc = Start.loc + Region.length
  Hap.ind = intersect(which(SNP.location>Start.loc),which(SNP.location<End.loc))
  return(Hap.ind)
}
