\name{Sample_Haplotypes}
\alias{Sample_Haplotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Choose random block of haplotypes for theoretical power}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sample_Haplotypes(SNP.location, Region.length, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNP.location}{
%%     ~~Describe \code{SNP.location} here~~
}
  \item{Region.length}{
%%     ~~Describe \code{Region.length} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (SNP.location, Region.length, n) 
{
    MAX.loc = max(SNP.location)
    MIN.loc = min(SNP.location)
    Start.loc = floor(runif(1) * (MAX.loc - MIN.loc - Region.length) + 
        MIN.loc)
    End.loc = Start.loc + Region.length
    Hap.ind = intersect(which(SNP.location > Start.loc), which(SNP.location < 
        End.loc))
    return(Hap.ind)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
