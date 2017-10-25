\name{CEPSKAT-package}
\alias{CEPSKAT-package}
\alias{CEPSKAT}
\docType{package}
\title{SNP-set (Sequence) Kernel Association Test for Continuous Extreme Phenotype Samples}
\description{
Test for an association between a SNP set and a continuous normally distributed phenotype in a sample selected to have only individuals with extreme phenotypes.}
\details{
\tabular{ll}{
Package: \tab CEPSKAT\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-04-26\cr
License: \tab GPL-3\cr
}
There are two primary functions:

1) SKAT_Null_Model_CEP(formula, c1, c2, delta = 0.001, gamma0 = "default", 
    MAXITERNUM = 20, data = NULL, n.Resampling = 0) 
	
2) CEP_Theoretical_Power(Haplotypes = NULL, SNP.location = NULL, Region.length = 3000, 
    nruns = 200, weights.beta = c(1, 25), n = 500, tails = 0.2, 
    MAF.cutoff = 0.03, prop.caus = 0.2, effect.size = 0.2, prop.pos = 1, 
    alpha.level = 10^(-6), optimal = FALSE)
	
Use CEP_Theoretical_Power to compute theoretical power over a wide variety of user specified settings given
a set of haplotypes and their basepair location.
Use SKAT_Null_Model_CEP in conjunction with the SKAT package to find p-values for associations between
a SNP-set and a continuous phenotype when sampling from phenotypic extremes.
}
\author{
Ian Barnett
Maintainer: ibarnett@hsph.harvard.edu
}
\references{
Wu, M. C.*, Lee, S.*, Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93. \ * contributed equally.

~ EPS paper citation
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
}
\examples{
# EXAMPLE input (unrealistic but only for example use)
N.avail = 10000
P = 1000
SNP.loc = 1:P
Haplotypes = matrix(rbinom(N.avail*P,1,.003),nrow=N.avail)

# Theoretical power calculation
CEP_Theoretical_Power(Haplotypes,SNP.loc,Region.length=50,nruns=100,prop.caus=0.1,prop.pos=0.8,n=500,weights.beta=c(1,25),tails=0.25,effect.size=0.3,alpha.level = 10^(-6),optimal=TRUE)


# Example SKAT_Null_Modek_CEP usage with no association
G = cbind(rbinom(1000,2,.2),rbinom(1000,2,.15))
Y = rnorm(1000) # No association
obj=SKAT_Null_Model_CEP(Y~1,0,0)
SKAT(G,obj)$p.value


# Example SKAT_Null_Modek_CEP usage with high association
Y = G%*%c(.4,.4)+rnorm(1000) # High association
obj=SKAT_Null_Model_CEP(Y~1,0,0)
SKAT(G,obj)$p.value
}
