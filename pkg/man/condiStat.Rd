\name{condiStat}
\alias{condiStat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
conditional statistics
}
\description{
 given a multivariate normal with mean mu and
 covariance matrix sig, compute the conditional
 mean and variance if we fix the variables indexed
 with indi to be equal to z

}
\usage{
condiStat(sig, mu, z, indi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sig}{
%%     ~~Describe \code{sig} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{indi}{
%%     ~~Describe \code{indi} here~~
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
function( sig, mu, z, indi )
{
    m  <- nrow( sig ); jndi <- ( 1:m )[-indi]
    mv <- rep( 0, m ); mv[indi] <- mu[indi] - z; mv[jndi] <- mu[jndi]
    aa <- sweeper( cbind( sig, mv ), indi )
    return( list( cmean = aa[jndi,m+1], cdisp = aa[jndi,jndi] ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
