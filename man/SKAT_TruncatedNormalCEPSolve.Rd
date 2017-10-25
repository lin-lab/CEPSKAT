\name{SKAT_TruncatedNormalCEPSolve}
\alias{SKAT_TruncatedNormalCEPSolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Solve score equations for truncated normal}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SKAT_TruncatedNormalCEPSolve(formula, c1, c2, delta = 0.001, gamma0 = "default", MAXITERNUM = 20, data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{c1}{
%%     ~~Describe \code{c1} here~~
}
  \item{c2}{
%%     ~~Describe \code{c2} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{gamma0}{
%%     ~~Describe \code{gamma0} here~~
}
  \item{MAXITERNUM}{
%%     ~~Describe \code{MAXITERNUM} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, c1, c2, delta = 0.001, gamma0 = "default", 
    MAXITERNUM = 20, data = NULL) 
{
    preerror = 0
    X <- model.matrix(formula, data = data)
    Y <- model.frame(formula, data = data)[, 1]
    if (is.character(gamma0) && gamma0 == "default") {
        out.lm <- lm(formula, data = data)
        sigma.ols = mean(abs(out.lm$residuals))
        c = (c1 - c2)/2
        gamma0 = c(out.lm$coefficients, sigma.ols)
    }
    n = nrow(X)
    K = ncol(X)
    iternum = 0
    sigma = gamma0[K + 1]
    alpha0 = gamma0[-(K + 1)]
    while (T) {
        Xa = (X \%*\% alpha0)[, 1]
        J = matrix(0, nrow = K, ncol = K)
        a1 <- (c1 - Xa)/sigma
        a2 <- (c2 - Xa)/sigma
        f1 <- dnorm(a1)
        f2 <- dnorm(a2)
        F1 <- pnorm(a1)
        F2 <- pnorm(a2)
        Y_Xa <- (Y - Xa)
        f12 <- f1 - f2
        F12 <- F2 + 1 - F1
        J[1:K, 1:K] <- t(X) \%*\% (X * (-1 + (a2 * f2 - a1 * f1)/F12 + 
            (f12/F12)^2))/sigma^2
        Jinv = solve(J)
        V = matrix(0, nrow = K, ncol = 1)
        V[1:K, 1] = colSums((X/sigma) * (Y_Xa/sigma - f12/F12))
        sigma_next = uniroot(function(sigmaT) sum(-1/sigmaT + 
            (Y_Xa)^2/sigmaT^3 + ((a2/sigmaT) * f2 - ((c1 - X \%*\% 
            alpha0)/sigmaT^2) * f1)/F12), c(sigma/20, 10 * sigma))$root
        alphaNext = alpha0 - Jinv \%*\% V
        iternum = iternum + 1
        curerror = sum(abs(alphaNext - alpha0)) + abs(sigma - 
            sigma_next)
        if (iternum > MAXITERNUM && curerror > preerror) {
            stop("Newton-Raphson diverging. Try different initial guess for gamma0.")
        }
        preerror = curerror
        if (curerror < delta) {
            end = TRUE
            if (end) {
                MLEs = c(alphaNext, sigma)
                MLEs = data.frame(MLEs)
                rownames = rep(0, K + 1)
                for (i in 1:K) {
                  rownames[i] = paste("alpha", i - 1, sep = "")
                }
                rownames[1] = "intercept"
                rownames[K + 1] = "sigma"
                row.names(MLEs) = rownames
                alpha0 = alphaNext
                Xa = (X \%*\% alpha0)[, 1]
                re <- SKAT_TruncatedNormalCEP_Get_Moment(c1, 
                  c2, sigma, Xa)
                return(list(coef = MLEs, M = re$M, V = re$V, 
                  W = re$W, Mu = re$Mu, Y = Y, X1 = X, sigma = sigma, 
                  Xa = Xa))
            }
        }
        alpha0 = alphaNext
        sigma = sigma_next
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
