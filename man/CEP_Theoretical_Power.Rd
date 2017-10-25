\name{CEP_Theoretical_Power}
\alias{CEP_Theoretical_Power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Theoretical Power Calculation for CEP-SKAT}
\description{Given haplotype data, this function will provide theoretical power calculations over a variety of user-specified settings.}
\usage{
CEP_Theoretical_Power(Haplotypes = NULL, SNP.location = NULL, Region.length = 3000, nruns = 200, weights.beta = c(1, 25), n = 500, tails = 0.2, MAF.cutoff = 0.03, prop.caus = 0.2, effect.size = 0.2, prop.pos = 1, alpha.level = 10^(-6), optimal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Haplotypes}{N.avail by P matrix, where N.avail=# of haplotypes and P=# of mutations.}
  \item{SNP.location}{Vector of base pair locations for Haplotypes (assumed to be on the same chromsome). Vector is of length P.}
  \item{Region.length}{Length (in base pairs) of region of interest.}
  \item{nruns}{Number of random regions selected. Estimated power is averaged over all nruns regions, so increased nruns improves the accuracy of the power estimate (but leads to longer computation)}
  \item{weights.beta}{The parameters of the beta distribution which will be used as weights in the K matrix. Weights will be generated using the Beta distribution function with these parameters (default is alpha=1, beta=25) evaluated at the MAF of the variant.}
  \item{n}{Sample size.}
  \item{tails}{This should be a number in the range (0,0.5) that specifies the proportion of the study population that will be selected from each tail of the phenotype distribution to be in the sample of size n.}
  \item{MAF.cutoff}{Variants with MAF > MAF.cutoff are not ever considered to be causal.}
  \item{prop.caus}{The proportion of variants that will be considered causal.}
  \item{effect.size}{This tuning parameter controls the effect size of all causal variants. The effect sizes are calculated as -effect.size*log10(MAF).}
  \item{prop.pos}{The proportion of causal variants that have a positive effect on phenotype. The other causal variants are assumed to have a negative effect.}
  \item{alpha.level}{The significance level of the test. Defaults to 10^-6 to mirror the scale of genome-wide association studies.}
  \item{optimal}{A TRUE value will use optimal SKAT. The default FALSE value will not use regular SKAT without the optimality adjustment.}
}
\details{
This function requires haplotype data as well as the SNP base pair location for each SNP. 
}
\value{This function returns a single estimate for power under the conditions given by the parameters.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Ian Barnett}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
N.avail = 10000
P = 1000
SNP.loc = 1:P
Haplotypes = matrix(rbinom(N.avail*P,1,.003),nrow=N.avail)
CEP_Theoretical_Power(Haplotypes,SNP.loc,Region.length=50,nruns=100,prop.caus=0.1,prop.pos=0.8,n=500,weights.beta=c(1,25),tails=0.25,effect.size=0.3,alpha.level = 10^(-6),optimal=TRUE)

## The function is currently defined as
function (Haplotypes = NULL, SNP.location = NULL, Region.length = 3000, 
    nruns = 200, weights.beta = c(1, 25), n = 500, tails = 0.2, 
    MAF.cutoff = 0.03, prop.caus = 0.2, effect.size = 0.2, prop.pos = 1, 
    alpha.level = 10^(-6), optimal = FALSE) 
{
    P = ncol(Haplotypes)
    N.avail = nrow(Haplotypes)
    N = as.integer(n/(2 * tails))
    MAF.full = colMeans(Haplotypes)
    if (length(which(MAF.full > 0.5)) != 0) {
        Haplotypes = Make_All_Minor_Alleles(Haplotypes, genotypes = FALSE, 
            MAF.full)
        MAF.full = colMeans(Haplotypes)
    }
    if (optimal) {
        rho = prop.caus^2 * (2 * prop.pos - 1)^2
    }
    power = rep(0, nruns)
    for (i in 1:nruns) {
        Hap.ind = Sample_Haplotypes(SNP.location, Region.length, 
            n)
        p = length(Hap.ind)
        H1 = sample(1:N.avail, N, replace = TRUE)
        H2 = sample(1:N.avail, N, replace = TRUE)
        G = as.matrix(Haplotypes[H1, Hap.ind] + Haplotypes[H2, 
            Hap.ind])
        MAF = colMeans(G)/2
        if (length(which(MAF > 0.5)) != 0) {
            G = Make_All_Minor_Alleles(G, genotype = TRUE, MAF)
            MAF = colMeans(G)/2
        }
        Beta = Get_Betas(MAF, MAF.cutoff, prop.caus, effect.size, 
            prop.pos)
        Mu = G \%*\% Beta
        Y = Mu + rnorm(nrow(G))
        cuts = quantile(Y, c(tails, 1 - tails))
        probSelect = pnorm(cuts[1], mean = Mu) + pnorm(cuts[2], 
            mean = Mu, lower.tail = FALSE)
        probSelect = probSelect * n/sum(probSelect)
        MAF = colSums(as.numeric(probSelect) * G)/sum(probSelect)
        W = Get_W_Beta_Dist(MAF, weights.beta)
        if (optimal) {
            W = W \%*\% (diag(1 - rho, nrow(W)) + matrix(rho, nrow = nrow(W), 
                ncol = nrow(W)))
        }
        G = as.numeric(sqrt(probSelect)) * G
        Mu = G \%*\% Beta
        tG.Mu = t(G) \%*\% Mu
        A = t(G) \%*\% G \%*\% W
        B = tG.Mu \%*\% t(tG.Mu) \%*\% W
        A.2 = A \%*\% A
        A.3 = A.2 \%*\% A
        A.4 = A.3 \%*\% A
        c.N = rep(0, 4)
        c.A = rep(0, 4)
        c.N[1] = sum(diag(A))
        c.N[2] = sum(diag(A.2))
        c.N[3] = sum(diag(A.3))
        c.N[4] = sum(diag(A.4))
        c.A[1] = sum(diag(B)) + c.N[1]
        c.A[2] = sum(diag(A \%*\% B)) + c.N[2]
        c.A[3] = sum(diag(A.2 \%*\% B)) + c.N[3]
        c.A[4] = sum(diag(A.3 \%*\% B)) + c.N[4]
        mu.Q = c.N[1]
        sigma.Q = sqrt(2 * c.N[2])
        if (is.na(sigma.Q)) {
            if (i == 1) {
                power[i] = 0
            }
            else {
                power[i] = power[i - 1]
            }
            next
        }
        s1 = c.N[3]/(c.N[2]^(3/2))
        s2 = c.N[4]/(c.N[2]^2)
        if (s1^2 > s2) {
            a = 1/(s1 - sqrt(s1^2 - s2))
            delta = s1 * (a^3) - a^2
        }
        else {
            a = 1/sqrt(s2)
            delta = 0
        }
        l = a^2 - 2 * delta
        mu.X = l + delta
        sigma.X = sqrt(2) * sqrt(l + 2 * delta)
        q.c = (qchisq(1 - alpha.level, df = l, ncp = delta) - 
            mu.X) * sigma.Q/sigma.X + mu.Q
        mu.Q = c.A[1]
        sigma.Q = sqrt(2 * c.A[2])
        if (is.na(sigma.Q)) {
            if (i == 1) {
                power[i] = 0
            }
            else {
                power[i] = power[i - 1]
            }
            next
        }
        s1 = c.A[3]/(c.A[2]^(3/2))
        s2 = c.A[4]/(c.A[2]^2)
        if (s1^2 > s2) {
            a = 1/(s1 - sqrt(s1^2 - s2))
            delta = s1 * (a^3) - a^2
        }
        else {
            a = 1/sqrt(s2)
            delta = 0
        }
        l = a^2 - 2 * delta
        mu.X = l + delta
        sigma.X = sqrt(2) * sqrt(l + 2 * delta)
        power[i] = pchisq(sigma.X * (q.c - mu.Q)/sigma.Q + mu.X, 
            df = max(l, 0), ncp = max(delta, 0), lower.tail = FALSE)
        if (floor(i/10) * 10 == i) {
            msg <- sprintf("\%d/\%d", i, nruns)
            print(msg)
        }
    }
    return(mean(power))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
