// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// poisson_disc_sampling
NumericMatrix poisson_disc_sampling(double radius, double size_x, double size_y, int num_samples_before_rejection);
RcppExport SEXP _poissonDiscSampling_poisson_disc_sampling(SEXP radiusSEXP, SEXP size_xSEXP, SEXP size_ySEXP, SEXP num_samples_before_rejectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< double >::type size_x(size_xSEXP);
    Rcpp::traits::input_parameter< double >::type size_y(size_ySEXP);
    Rcpp::traits::input_parameter< int >::type num_samples_before_rejection(num_samples_before_rejectionSEXP);
    rcpp_result_gen = Rcpp::wrap(poisson_disc_sampling(radius, size_x, size_y, num_samples_before_rejection));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_poissonDiscSampling_poisson_disc_sampling", (DL_FUNC) &_poissonDiscSampling_poisson_disc_sampling, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_poissonDiscSampling(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
