// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gaussian_ridge
List gaussian_ridge(const Eigen::VectorXf& y, const Eigen::MatrixXf& X, double lambda, int iter_max, double tolerance);
RcppExport SEXP _ridges_gaussian_ridge(SEXP ySEXP, SEXP XSEXP, SEXP lambdaSEXP, SEXP iter_maxSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXf& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXf& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian_ridge(y, X, lambda, iter_max, tolerance));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ridges_gaussian_ridge", (DL_FUNC) &_ridges_gaussian_ridge, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_ridges(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}