\name{SKAT_TruncatedNormalCEP_Get_Moment}
\alias{SKAT_TruncatedNormalCEP_Get_Moment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get moments of Truncated Normal Distribution}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SKAT_TruncatedNormalCEP_Get_Moment(c1, c2, sigma, Xa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c1}{
%%     ~~Describe \code{c1} here~~
}
  \item{c2}{
%%     ~~Describe \code{c2} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{Xa}{
%%     ~~Describe \code{Xa} here~~
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
function (c1, c2, sigma, Xa) 
{
    a1 <- (c1 - Xa)/sigma
    a2 <- (c2 - Xa)/sigma
    f1 <- dnorm(a1)
    f2 <- dnorm(a2)
    F1 <- pnorm(a1)
    F2 <- pnorm(a2)
    f12 <- f2 - f1
    F12 <- F2 + 1 - F1
    M <- f12/F12 * sigma
    V <- (a2 * f2 - a1 * f1)/F12 + M^2/sigma^2
    W <- (1 - V) * sigma^2
    Mu <- Xa - M
    return(list(M = M, V = V, W = W, Mu = Mu))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
