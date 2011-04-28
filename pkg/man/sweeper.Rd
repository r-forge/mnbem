\name{sweeper}
\alias{sweeper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beaton's sweep function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sweeper(a, indi, rev = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{indi}{
%%     ~~Describe \code{indi} here~~
}
  \item{rev}{
%%     ~~Describe \code{rev} here~~
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
function( a, indi, rev=FALSE ){
    d = dim( a )
    if( d[1]!=d[2] ) stop( "Only works on square matrix" )
    if( rev ){ rv = -1 } else { rv = 1 }
    z <- .Call( "SWEEP",
                as.integer( indi-1 ),
                as.double( a ),
                as.integer( d[1] ),
                as.double( rv )
              );
    return( array( z, d ) );
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
