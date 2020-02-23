// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// encode_alphabet
Rcpp::NumericVector encode_alphabet(const Rcpp::NumericVector& input);
RcppExport SEXP _seqR_encode_alphabet(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(encode_alphabet(input));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_seqR_encode_alphabet", (DL_FUNC) &_seqR_encode_alphabet, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_seqR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
