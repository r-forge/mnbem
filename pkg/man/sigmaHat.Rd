\name{sigmaHat}
\alias{sigmaHat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate sigma hat
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigmaHat(vmat, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vmat}{
%%     ~~Describe \code{vmat} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
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
function( vmat, v )
{
    n <- nrow( v ); m <- nrow( vmat ) / n; mm <- 1:m; shat <- matrix( 0, m, m )
    for ( i in 1:n ){
         for ( k in 1:n ){
            shat <- shat + v[i,k] * vmat[mm + ( i - 1 ) * m, mm + ( k - 1 ) * m]
         }
    }
    return(shat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
