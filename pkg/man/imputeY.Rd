\name{imputeY}
\alias{imputeY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
impute Y.
}
\description{
 given a multivariate normal with mean mu and
 covariance matrix sig, and a vector of observations y,
 impute the missing values in y and the corresponding
 conditional covariance matrix

}
\usage{
imputeY(y, sig, mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{sig}{
%%     ~~Describe \code{sig} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
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
function( y, sig, mu )
{
    n <- length( y ); nm < -is.na( y )
    indi <- which( !nm ); v <- matrix( 0, n, n )
    cs <- condiStat( sig, mu, y[indi], indi )
    y[nm] <- cs$cmean; v[nm,nm] <- cs$cdisp
    return( list( yimp = y, vimp = v ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
