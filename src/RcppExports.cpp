// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>

using namespace Rcpp;

// rcpp_bin_quantile_rec
int rcpp_bin_quantile_rec(SEXP rn, SEXP rp, SEXP rnmin, SEXP rnmax, SEXP rminprob, SEXP rdir);
RcppExport SEXP motifDiverge_rcpp_bin_quantile_rec(SEXP rnSEXP, SEXP rpSEXP, SEXP rnminSEXP, SEXP rnmaxSEXP, SEXP rminprobSEXP, SEXP rdirSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type rn(rnSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rp(rpSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rnmin(rnminSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rnmax(rnmaxSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rminprob(rminprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rdir(rdirSEXP );
        int __result = rcpp_bin_quantile_rec(rn, rp, rnmin, rnmax, rminprob, rdir);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_convolve
arma::vec rcpp_convolve(SEXP rx, SEXP ry);
RcppExport SEXP motifDiverge_rcpp_convolve(SEXP rxSEXP, SEXP rySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type rx(rxSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ry(rySEXP );
        arma::vec __result = rcpp_convolve(rx, ry);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_prob_sym
NumericVector rcpp_prob_sym(SEXP rk, SEXP rn, SEXP rp10, SEXP rp01, SEXP rerr);
RcppExport SEXP motifDiverge_rcpp_prob_sym(SEXP rkSEXP, SEXP rnSEXP, SEXP rp10SEXP, SEXP rp01SEXP, SEXP rerrSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type rk(rkSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rn(rnSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rp10(rp10SEXP );
        Rcpp::traits::input_parameter< SEXP >::type rp01(rp01SEXP );
        Rcpp::traits::input_parameter< SEXP >::type rerr(rerrSEXP );
        NumericVector __result = rcpp_prob_sym(rk, rn, rp10, rp01, rerr);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
