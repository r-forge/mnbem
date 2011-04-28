#include <R.h>
#include <Rinternals.h>

// sweep operator function
SEXP SWEEP( SEXP, SEXP, SEXP, SEXP );

// sweep operator
// indi: index to sweep, could be scaler or vector
// mat:  d by d square matrix to sweep needs to be numeric, treated as double
// dim:  dimension of mat
// rev:  1/-1 depending on forward or reverse sweep

SEXP
SWEEP( SEXP indi, SEXP mat, SEXP dim, SEXP rev )
{
    int  j, p = length( indi ), n = length( mat ), d = INTEGER( dim )[0], nProtected = 0;
    double pv;
    SEXP SWAPMAT;
    PROTECT( SWAPMAT = allocVector( REALSXP, n * p ) ); ++nProtected;
    for( int k = 0; k < p; k++ ) {
        j = INTEGER( indi )[k];
        pv = REAL( mat )[( j * d ) + j];
        REAL( SWAPMAT )[( j * d ) + j ] = -1/pv;
        for( int r = 0; r < d; r++ ) {
            if( j != r ){
                REAL( SWAPMAT )[( r * d ) + j ] = REAL(rev)[0] * REAL( mat )[( r * d ) + j ]/pv;
            }
        }
        for( int c = 0; c < d; c++ ) {
            if( j != c ){
                REAL( SWAPMAT )[( j * d ) + c ] = REAL(rev)[0] * REAL( mat )[( j * d ) + c ]/pv;
            }
        }
        for( int r = 0; r < d; r++ ) {
            for( int c = 0; c < d; c++ ) {
                if( j != c && j != r ){
                    REAL( SWAPMAT )[( r * d ) + c ] = REAL( mat )[( r * d ) + c ]-REAL( mat )[( j * d ) + c] * REAL( mat )[( r * d ) + j]/pv;
                    //REAL( SWAPMAT )[( r * d ) + c ] =REAL( SWAPMAT )[( c * d ) + r ] ;
                }
            }
        }
    }
    UNPROTECT( nProtected );
    return SWAPMAT;
}

