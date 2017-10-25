CEP_Theoretical_Power <-
function(Haplotypes=NULL,SNP.location=NULL,Region.length=3000,nruns=200,weights.beta=c(1,25),n=500,tails=0.20,MAF.cutoff=0.03,prop.caus=0.2,effect.size=.2,prop.pos=1,alpha.level = 10^(-6),optimal=FALSE){
  P = ncol(Haplotypes)
  N.avail = nrow(Haplotypes)
  N = as.integer(n/(2*tails))
  MAF.full = colMeans(Haplotypes)
  if(length(which(MAF.full>0.5))!=0){
    Haplotypes = Make_All_Minor_Alleles(Haplotypes,genotypes=FALSE,MAF.full)
    MAF.full = colMeans(Haplotypes)
  }
  if(optimal){
    rho = prop.caus^2*(2*prop.pos-1)^2
  }
  power = rep(0,nruns)
  for(i in 1:nruns){
    Hap.ind = Sample_Haplotypes(SNP.location,Region.length,n)
    p = length(Hap.ind)
    H1 = sample(1:N.avail,N,replace=TRUE)
    H2 = sample(1:N.avail,N,replace=TRUE)
    G = as.matrix(Haplotypes[H1,Hap.ind]+Haplotypes[H2,Hap.ind])
    MAF = colMeans(G)/2
    if(length(which(MAF>0.5))!=0){
      G = Make_All_Minor_Alleles(G,genotype=TRUE,MAF)
      MAF = colMeans(G)/2
    }
    Beta = Get_Betas(MAF,MAF.cutoff,prop.caus,effect.size,prop.pos)
    Mu = G%*%Beta
    Y = Mu + rnorm(nrow(G))
    cuts = quantile(Y,c(tails,1-tails))
    probSelect = pnorm(cuts[1],mean=Mu)+pnorm(cuts[2],mean=Mu,lower.tail=FALSE)
    probSelect = probSelect*n/sum(probSelect)
    MAF = colSums(as.numeric(probSelect)*G)/sum(probSelect)             
    W = Get_W_Beta_Dist(MAF,weights.beta)
    if(optimal){
      W = W %*% (diag(1-rho,nrow(W))+matrix(rho,nrow=nrow(W),ncol=nrow(W)))
    }
    G = as.numeric(sqrt(probSelect))*G
    Mu = G%*%Beta
    tG.Mu = t(G)%*%Mu
    A = t(G)%*%G%*%W
    B = tG.Mu%*%t(tG.Mu)%*%W
    A.2 = A %*% A
    A.3 = A.2 %*% A
    A.4 = A.3 %*% A
    # Obtain cumulants under null (N) and alternative (A)
    c.N = rep(0,4)
    c.A = rep(0,4)
    c.N[1]= sum(diag(A))
    c.N[2]= sum(diag(A.2))
    c.N[3]= sum(diag(A.3))
    c.N[4]= sum(diag(A.4))
    c.A[1]= sum(diag(B))+c.N[1]
    c.A[2]= sum(diag(A %*% B))+c.N[2]
    c.A[3]= sum(diag(A.2 %*% B))+c.N[3]
    c.A[4]= sum(diag(A.3 %*% B))+c.N[4]
    # Obtain quantile under Null
    mu.Q = c.N[1]
    sigma.Q = sqrt(2*c.N[2])
    if(is.na(sigma.Q)){ # In case of an error, skips this iteration
      if(i==1){power[i]=0}
      else{power[i]=power[i-1]}
      next
    }
    s1 = c.N[3]/(c.N[2]^(3/2))
    s2 = c.N[4]/(c.N[2]^2)
    if(s1^2>s2){
      a = 1/(s1-sqrt(s1^2-s2))
      delta = s1*(a^3)-a^2
    } else{
      a = 1/sqrt(s2)
      delta = 0
    }
    l = a^2-2*delta
    mu.X = l+delta
    sigma.X = sqrt(2)*sqrt(l+2*delta)
    q.c = (qchisq(1-alpha.level,df=l,ncp=delta)-mu.X)*sigma.Q/sigma.X+mu.Q
    # Obtain power under Alternative
    mu.Q = c.A[1]
    sigma.Q = sqrt(2*c.A[2])
    if(is.na(sigma.Q)){ # In case of an error, skips this iteration
      if(i==1){power[i]=0}
      else{power[i]=power[i-1]}
      next
    }
    s1 = c.A[3]/(c.A[2]^(3/2))
    s2 = c.A[4]/(c.A[2]^2)
    if(s1^2>s2){
      a = 1/(s1-sqrt(s1^2-s2))
      delta = s1*(a^3)-a^2
    } else{
      a = 1/sqrt(s2)
      delta = 0
    }
    l = a^2-2*delta
    mu.X = l+delta
    sigma.X = sqrt(2)*sqrt(l+2*delta)
    power[i]=pchisq(sigma.X*(q.c-mu.Q)/sigma.Q+mu.X,df=max(l,0),ncp=max(delta,0),lower.tail=FALSE)
    if(floor(i/10) * 10 == i){
      msg<-sprintf("%d/%d",i,nruns)
      print(msg)
    }
  }
  return(mean(power))
}
