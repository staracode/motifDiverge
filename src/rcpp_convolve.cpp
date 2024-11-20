
//===========================
// CONVOLUTION OF TWO VECTORS 
//===========================

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec rcpp_convolve(SEXP rx, SEXP ry) {


	arma::vec x = Rcpp::as<arma::vec>(rx) ; // vector x to be convoluted
	arma::vec y = Rcpp::as<arma::vec>(ry) ; // vector y to be convoluted by
	arma::vec z = conv(x,y)               ;

	int nx = x.n_elem ;

	return(z);

}
