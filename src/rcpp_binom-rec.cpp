
#include <RcppGSL.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <cmath>

using namespace Rcpp;

//=========================================
// RECURSEIVE NESTED INTERVALS FOR BINOMIAL
//=========================================

int bin_quantile_rec(	int     n,
                        double  p,
                        double  eps,
                        int 	min,
                        int 	max,
                        int 	dir){

        if(max == min){
                return max;
        }
        if(max - min == 1){
                if(dir==1) return(max); //- right-side
		return(min);
        }
        else{
                int mid = int (ceil ( min + (max-min)/2.) ); // ceiling
                double val = gsl_ran_binomial_pdf((unsigned int) mid, p, (unsigned int) n);
                if(dir==1){
                        if(val > eps) return bin_quantile_rec(n,p,eps,mid,max,1);
                        else return bin_quantile_rec(n,p,eps,min,mid,1);
                }
                else{
                        if(val > eps) return bin_quantile_rec(n,p,eps,min,mid,-1);
                        else return bin_quantile_rec(n,p,eps,mid,max,-1);
                }
        }
}	

// [[Rcpp::export]]
int rcpp_bin_quantile_rec( SEXP rn, SEXP rp, SEXP rnmin, SEXP rnmax, SEXP rminprob, SEXP rdir) { 


	int n 		= as<int>(rn)		;
	double p  	= as<double>(rp)	;
	int nmin 	= as<int>(rnmin)	;
	int nmax 	= as<int>(rnmax)	;
	double minprob  = as<double>(rminprob)	;
	int dir 	= as<int>(rdir)		;

	int res		= dir > 0 ? nmax : nmin ;
	int mode        = int( floor((n+1)*p) > floor((n+1)*p-1) ? floor((n+1)*p) : floor((n+1)*p-1) ) ; 
	
	nmin = dir > 0 ? mode : nmin ;
	nmax = dir > 0 ? nmax : mode ;
	
	res = bin_quantile_rec(n,p,minprob,nmin,nmax,dir) ;

	return(res);

}


