SKAT_TruncatedNormalCEPSolve <-
function(formula, c1, c2, delta=0.001, gamma0="default", MAXITERNUM=20, data=NULL){
  preerror = 0
  X<-model.matrix(formula,data=data)
  Y<-model.frame(formula, data=data)[,1]
  if(is.character(gamma0) && gamma0=="default"){
    out.lm<-lm(formula, data=data)
    sigma.ols =mean(abs(out.lm$residuals))
    c = (c1-c2)/2 
    gamma0 = c(out.lm$coefficients,sigma.ols)
  }
  n = nrow(X)
  K = ncol(X)
  iternum = 0
  sigma = gamma0[K+1]
  alpha0 = gamma0[-(K+1)]
  while(T){
    Xa = (X%*%alpha0)[,1]
    J = matrix(0,nrow=K,ncol=K)
    a1<-(c1-Xa)/sigma
    a2<-(c2-Xa)/sigma
    f1<-dnorm(a1)
    f2<-dnorm(a2)
    F1<-pnorm(a1)
    F2<-pnorm(a2)
    Y_Xa<-(Y-Xa)
    f12<-f1-f2
    F12<-F2 + 1 - F1
    J[1:K, 1:K]<-t(X) %*% ( X * (-1+(a2*f2-a1*f1)/F12+(f12/F12)^2)) / sigma^2
    Jinv = solve(J)
    V = matrix(0,nrow=K,ncol=1)
    V[1:K,1] = colSums((X/sigma)*(Y_Xa/sigma-f12/F12))
    sigma_next = uniroot(function(sigmaT) sum(-1/sigmaT+(Y_Xa)^2/sigmaT^3+((a2/sigmaT)*f2-((c1-X%*%alpha0)/sigmaT^2)*f1)/F12),c(sigma/20,10*sigma))$root
    alphaNext = alpha0-Jinv%*%V
    iternum = iternum + 1
    curerror = sum(abs(alphaNext-alpha0))+abs(sigma-sigma_next)
    if(iternum > MAXITERNUM && curerror > preerror){stop("Newton-Raphson diverging. Try different initial guess for gamma0.")}
    preerror = curerror  
    if(curerror < delta){
      end = TRUE
      if(end){
        MLEs = c(alphaNext,sigma)
        MLEs = data.frame(MLEs)
        rownames = rep(0,K+1)
        for(i in 1:K){
          rownames[i] = paste("alpha",i-1,sep="")
        }
        rownames[1] = "intercept"
        rownames[K+1] = "sigma"
        row.names(MLEs) = rownames
        alpha0 = alphaNext
        Xa = (X%*%alpha0)[,1]
        
        re<-SKAT_TruncatedNormalCEP_Get_Moment(c1, c2, sigma, Xa)
        return(list(coef=MLEs, M = re$M, V= re$V, W=re$W, Mu=re$Mu, Y=Y, X1=X, sigma=sigma, Xa=Xa))
      }
    }
    alpha0 = alphaNext
    sigma=sigma_next
  }
}
