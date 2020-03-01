#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>

#include "alphabet_encoder.h"
#include <memory>
#include <string>

typedef short ENCODED_ELEM_T;

auto getEncoding(Rcpp::IntegerVector& input) {
  return getAlphabetEncoding<Rcpp::IntegerVector, int, int, ENCODED_ELEM_T>(
    input,
    [](const int& elem) -> int { return elem; }
  );
}
  
auto getEncoding(Rcpp::NumericVector& input) {
  return getAlphabetEncoding<Rcpp::NumericVector, double, double, ENCODED_ELEM_T>(
    input,
    [](const double& elem) -> double { return elem; }
  );
}

auto getEncoding(Rcpp::StringVector& input) {
  return getAlphabetEncoding<Rcpp::StringVector, Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T>(
    input,
    [](const Rcpp::String::StringProxy& elem) -> std::string { return Rcpp::as<std::string>(elem); }
  );
}

#endif //CONCRETE_ENCODERS_H
