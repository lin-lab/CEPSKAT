\name{SKAT_Null_Model_CEP}
\alias{SKAT_Null_Model_CEP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit the Null model for CEP-SKAT}
\description{
This function should be used in conjunction with the SKAT function of the SKAT package. A null model is fit and
prepared to be passed to the SKAT function in order to perform the score test for association.
}
\usage{
SKAT_Null_Model_CEP(formula, c1, c2, delta = 0.001, gamma0 = "default", MAXITERNUM = 20, data = NULL, n.Resampling = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{The null model to be fit.}
  \item{c1}{This is the smallest (lower bound) for the high extreme phenotypes. In other words c1 is the lower bound for the phenotypes sampled from the upper tail of the distribution.}
  \item{c2}{This is the largest (upper bound) for the low extreme phenotypes. In other words c2 is the upper bound for the phenotypes sampled from the lower tail of the distribution.}
  \item{delta}{Specifies speed and accuracy of convergence. Large values of delta converge faster but are less accurate, small values of delta converge slower but are more accurate.}
  \item{gamma0}{The initial guess for the MLEs. Default uses estimates from OLS.}
  \item{MAXITERNUM}{If the algorithm has not converged after this many iterations, the function exits. This protects against divergence of the Newton-Raphson algorithm.}
  \item{data}{Original data.frame with phenotype and covariates. Optional: only use if formula contains references to items in data.}
  \item{n.Resampling}{Defaults to 0. If resampling is to be used, this specifies how many times to resample.}
}
\details{
SKAT tests for an association between a SNP-set and phenotype through a score test where the model is
fit under the null hypothesis of no association. This function fits this null model that is independent of
genotype data. After fitting this model, the SKAT package can handle the output combined with the genotypes
to produce p-values.
}
\value{
An item of class "SKAT_Null_Model" is returned which can be interpreted by the SKAT function in the SKAT package.}
\references{
Wu, M. C.*, Lee, S.*, Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93. \ * contributed equally.
}
\author{
Ian Barnett
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example with no association
G = cbind(rbinom(1000,2,.2),rbinom(1000,2,.15))
Y = rnorm(1000) # No association
ind_sample = order(Y)[c(1:250,751:1000)]
Y_extreme = Y[ind_sample]
G_extreme = G[ind_sample,]
c1 = min(Y[ind_sample[-(1:250)]])
c2 = max(Y[ind_sample[1:250]])
obj=SKAT_Null_Model_CEP(Y_extreme~1,c1,c2)
SKAT(G_extreme,obj)$p.value

# Example with high association
Y = G%*%c(.4,.4)+rnorm(1000) # High association
ind_sample = order(Y)[c(1:250,751:1000)]
Y_extreme = Y[ind_sample]
G_extreme = G[ind_sample,]
c1 = min(Y[ind_sample[-(1:250)]])
c2 = max(Y[ind_sample[1:250]])
obj=SKAT_Null_Model_CEP(Y_extreme~1,c1,c2)
SKAT(G_extreme,obj)$p.value



## The function is currently defined as
function (formula, c1, c2, delta = 0.001, gamma0 = "default", 
    MAXITERNUM = 20, data = NULL, n.Resampling = 0) 
{
    out <- SKAT_TruncatedNormalCEPSolve(formula, c1, c2, delta, 
        gamma0, MAXITERNUM, data)
    res = out$Y - out$Mu
    pi_1 = out$W
    X1 = out$X1
    n <- length(res)
    id_include = 1:n
    res.out = NULL
    type.Resampling = NULL
    if (n.Resampling > 0) {
        res.out <- SKAT_Get_TruncatedNormal(length(res), out$Xa, 
            out$sigma, c1, c2, nSet = n.Resampling) - out$Mu
        res.out <- t(t(res.out) - colMeans(res.out))
    }
    re <- list(res = res, X1 = X1, res.out = res.out, out_type = "D", 
        n.Resampling = n.Resampling, type.Resampling = NULL, 
        id_include = id_include, mu = out$Mu, pi_1 = pi_1, coef = out$coef, 
        Xa = out$Xa)
    class(re) <- "SKAT_NULL_Model"
    return(re)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
