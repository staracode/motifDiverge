
#include <RcppGSL.h>
#include <gsl/gsl_sf_gamma.h>
#include <cmath>

using namespace Rcpp;

//================
// CALCULATE P-SYM
//================

// [[Rcpp::export]]
NumericVector rcpp_prob_sym( SEXP rk, SEXP rn, SEXP rp10, SEXP rp01, SEXP rerr) { 


	int 	 n = as<int>(rn)	;
	int	 k = as<int>(rk)	;
	double p10 = as<double>(rp10)	;
	double p01 = as<double>(rp01)	;
	double err = as<double>(rerr)	;

	NumericVector resVec(3); 

	//- positive vs. negative k case
	//------------------------------
	if(k<0){
        	k = abs(k);
        	//- swap p10 and p01
        	p10 = p10+p01;
        	p01 = p10-p01;
        	p10 = p10-p01;
	}


	int    imax   = (int) floor( (n-k)/2.);
	double cerr   = INFINITY; // current error
	double lcnst  = gsl_sf_lnfact((unsigned int) n) +  k*log(p10) + (n-k)*log(1-p10-p01);
	double lsmnd1 = (-1)*gsl_sf_lnfact((unsigned int) k) - gsl_sf_lnfact((unsigned int) (n-k));
	double lsmnd2 = 0.;
	double res    = exp(lcnst+lsmnd1);
	int    ntrms  = 0;

	for(int i=1; i<=imax; i++){
		// increment
        	lsmnd2  = lsmnd1 + (( log(p01) + log(p10) - 2*log(1-p10-p01) )) - (( log(k+i) + log(i) - log(n-k-2*(i-1)) - log(n-k-2*(i-1)-1) )) ;
        	res    += exp(lcnst + lsmnd2);
 		// caclulate error bound
        	if(lsmnd2 < lsmnd1){
               		if( i < imax){
                       		cerr = exp(lcnst+lsmnd1) * ( (1-exp((lsmnd2-lsmnd1)*(imax-i)))/(1-exp(lsmnd2-lsmnd1)) -1);
               		} else {
                       		// summed up all terms
                       		cerr = 0;
		       		ntrms = imax;
               		}
        	}	
		// are we done?
       		if(cerr < err){
			ntrms = i; 
			break;
       		}
       		// next summand
       		lsmnd1 = lsmnd2;

	}

	resVec[0] = res;
	resVec[1] = cerr;
	resVec[2] = ntrms; 

	return(resVec);

}


