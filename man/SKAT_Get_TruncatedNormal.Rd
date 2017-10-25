\name{SKAT_Get_TruncatedNormal}
\alias{SKAT_Get_TruncatedNormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get Truncated Normal MLEs}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SKAT_Get_TruncatedNormal(n, mu, sigma, c1, c2, nSet, nUpper = n/2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{c1}{
%%     ~~Describe \code{c1} here~~
}
  \item{c2}{
%%     ~~Describe \code{c2} here~~
}
  \item{nSet}{
%%     ~~Describe \code{nSet} here~~
}
  \item{nUpper}{
%%     ~~Describe \code{nUpper} here~~
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
function (n, mu, sigma, c1, c2, nSet, nUpper = n/2) 
{
    prob1 <- 1 - pnorm(c1, mean = mu, sd = sigma)
    prob2 <- pnorm(c2, mean = mu, sd = sigma)
    prob <- prob1/(prob1 + prob2)
    temp.b <- rbinom(n * nSet, 1, rep(prob, nSet))
    lower = rep(-10^5, n * nSet) * (1 - temp.b) + c1
    upper = rep(10^5, n * nSet) * temp.b + c2
    re <- rtnorm(n * nSet, mean = mu, sd = sigma, lower = lower, 
        upper = upper)
    re.M <- matrix(re, ncol = nSet, byrow = FALSE)
    cutoff <<- cbind(lower, upper, temp.b, mu)
    rm(re)
    rm(lower)
    rm(upper)
    return(re.M)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
