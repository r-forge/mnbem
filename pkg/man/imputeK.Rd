\name{imputeK}
\alias{imputeK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
impute K
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imputeK(y, v, w, mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
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
function( y, v, w, mu )
{
    n <- nrow( y ); m <- ncol( y )
    vsig<- kronecker( v, w ); vy <- as.vector( t( y ) ); vmu <- as.vector( t( mu ) )
    vim <- imputeY( vy, vsig, vmu )
    return( list( vim = matrix( vim$yimp, n, m, byrow = TRUE ), vvv = vim$vimp ) )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
