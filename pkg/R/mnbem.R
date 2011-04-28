#
#   mnem package
#   Copyright (C) 2009  Jan de Leeuw <deleeuw@stat.ucla.edu>
#   UCLA Department of Statistics, Box 951554, Los Angeles, CA 90095-1554
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warrant(y) of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################
#
# version 0.1, 2009-06-06   Initial  Release
# version 0.2, 2009-06-08   Added some stuff
#

# which(outer(i,j,function(x,y) is.na(x)&is.na(y)))

imputeK <- function( y, v, w, mu )
{
    n <- nrow( y ); m <- ncol( y )
    vsig<- kronecker( v, w ); vy <- as.vector( t( y ) ); vmu <- as.vector( t( mu ) )
    vim <- imputeY( vy, vsig, vmu )
    return( list( vim = matrix( vim$yimp, n, m, byrow = TRUE ), vvv = vim$vimp ) )
}

sigmaHat <- function( vmat, v )
{
    n <- nrow( v ); m <- nrow( vmat ) / n; mm <- 1:m; shat <- matrix( 0, m, m )
    for ( i in 1:n ){
         for ( k in 1:n ){
            shat <- shat + v[i,k] * vmat[mm + ( i - 1 ) * m, mm + ( k - 1 ) * m]
         }
    }
    return(shat)
}

gammaHat <- function( vmat, w )
{
    m <- nrow( w ); n <- nrow( vmat ) / m; mm <- ( 0:( n - 1 ) ) * m; ghat <- matrix( 0, n, n )
    for( j in 1:m ) {
        for( l in 1:m ){ 
            ghat <- ghat + w[ j, l ] * vmat[ mm + j, mm + l ]
        }
    }
return(ghat)
}

# given a multivariate normal with mean mu and
# covariance matrix sig, and a vector of observations y,
# impute the missing values in y and the corresponding
# conditional covariance matrix

imputeY <- function( y, sig, mu )
{
    n <- length( y ); nm <- is.na( y )
    indi <- which( !nm ); v <- matrix( 0, n, n )
    cs <- condiStat( sig, mu, y[indi], indi )
    y[nm] <- cs$cmean; v[nm,nm] <- cs$cdisp
    return( list( yimp = y, vimp = v ) )
}

# given a multivariate normal with mean mu and
# covariance matrix sig, compute the conditional
# mean and variance if we fix the variables indexed
# with indi to be equal to z

condiStat <- function( sig, mu, z, indi )
{
    m  <- nrow( sig ); jndi <- ( 1:m )[-indi]
    mv <- rep( 0, m ); mv[indi] <- mu[indi] - z; mv[jndi] <- mu[jndi]
    aa <- sweeper( cbind( sig, mv ), indi )
    return( list( cmean = aa[jndi,m+1], cdisp = aa[jndi,jndi] ) )
}

# beaton's sweep function
#dyn.load("mnbem.so")
sweeper <-function( a, indi, rev=FALSE ){
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
#
#sweeper<-function(a,indi) {
#n<-nrow(a); m<-length(indi)
#for (j in indi) {
#    pv<-a[j,j]
#    if (pv == 0) next()
#    pr<-a[j,-j]; pc<-a[-j,j]
#    a[j,j] <- -1/pv
#    a[j,-j] <- pr/pv
#    a[-j,j] <- pc/pv
#    a[-j,-j] <- a[-j,-j]-outer(pc,pr)/pv
#    }
#return(a)
#}

# this minimizes trace{(y-mu)'(y-mu)} over the
# missing elements of y 

imputeOLS <- function( y, mu ){
     return( ifelse( is.na( y ), mu, y ) )
}

imputeWLS<-function(y,k,ginv,sinv,mu,eps=1e-6,itmax=100,verbose=FALSE)
{
n<-nrow(y); m<-ncol(y); itel<-1
res<-y-mu; ras<-ginv%*%res%*%sinv; ff<-sum(ras*res)
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
