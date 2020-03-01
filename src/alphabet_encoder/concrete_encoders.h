#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>

#include "alphabet_encoder.h"
#include <memory>

typedef int ENCODED_T;

auto getEncoding(const Rcpp::IntegerVector& input) {
  return getAlphabetEncoding<Rcpp::IntegerVector, int, int, ENCODED_T>(
    input,
    [](const int& elem) -> int { return elem; }
  );
}

auto getEncoding(const Rcpp::NumericVector& input) {
  return getAlphabetEncoding<Rcpp::NumericVector, double, double, ENCODED_T>(
    input,
    [](const double& elem) -> double { return elem; }
  );
}

#endif //CONCRETE_ENCODERS_H