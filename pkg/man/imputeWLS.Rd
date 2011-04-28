\name{imputeWLS}
\alias{imputeWLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
impute using WLS
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imputeWLS(y, k, ginv, sinv, mu, eps = 1e-06, itmax = 100, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{ginv}{
%%     ~~Describe \code{ginv} here~~
}
  \item{sinv}{
%%     ~~Describe \code{sinv} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{itmax}{
%%     ~~Describe \code{itmax} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function(y,k,ginv,sinv,mu,eps=1e-6,itmax=100,verbose=FALSE)
{
n<-nrow(y); m<-ncol(y); itel<-1
res<-y-mu; ras<-ginv\%*\%res\%*\%sinv; ff<-sum(ras*res)
wgt<-outer(diag(ginv),diag(sinv))
repeat {
    thmax<-0
    for (i in 1:n) for (j in 1:m) {
	    if (!is.na(k[i,j])) next()
	    rij<-ras[i,j]; wij<-wgt[i,j]
        th<--rij/wij
        thmax<-max(thmax,abs(th))
        z[i,j]<-z[i,j]+th
        ff<-ff-(rij^2)/wij
        ras<-ras+th*outer(ginv[,i],sinv[,j])
        }
    if (verbose)
           cat("itel",formatC(itel,format="d",width=4)," maxth",formatC(thmax,format="f",digits=8,width=15)," func",formatC(ff,format="f",digits=8,width=15),"\n")
    if ((thmax < eps) || (itel == itmax)) break()
    itel<-itel+1
    }
return(y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
