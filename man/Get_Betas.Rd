\name{Get_Betas}
\alias{Get_Betas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get effect sizes for causal variants}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Get_Betas(MAF, MAF.cutoff, prop.caus, c, prop.pos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAF}{
%%     ~~Describe \code{MAF} here~~
}
  \item{MAF.cutoff}{
%%     ~~Describe \code{MAF.cutoff} here~~
}
  \item{prop.caus}{
%%     ~~Describe \code{prop.caus} here~~
}
  \item{c}{
%%     ~~Describe \code{c} here~~
}
  \item{prop.pos}{
%%     ~~Describe \code{prop.pos} here~~
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
function (MAF, MAF.cutoff, prop.caus, c, prop.pos) 
{
    p = length(MAF)
    Beta = rep(0, p)
    for (i in 1:p) {
        if (MAF[i] < MAF.cutoff && rnorm(1) < prop.caus && MAF[i] > 
            0) {
            Beta[i] = -c * log10(MAF[i]) * (2 * rbinom(1, 1, 
                prop.pos) - 1)
        }
    }
    return(Beta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
