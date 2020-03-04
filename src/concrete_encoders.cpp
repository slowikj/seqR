// [[Rcpp::plugins("c++17")]]
#include "concrete_encoders.h"

AlphabetEncoding<int, int, ENCODED_ELEM_T> getEncoding(Rcpp::IntegerVector& input) {
  return getAlphabetEncoding<Rcpp::IntegerVector, int, int, ENCODED_ELEM_T>(
      input,
      [](const int& elem) -> int { return elem; }
  );
}

AlphabetEncoding<double, double, ENCODED_ELEM_T> getEncoding(Rcpp::NumericVector& input) {
  return getAlphabetEncoding<Rcpp::NumericVector, double, double, ENCODED_ELEM_T>(
      input,
      [](const double& elem) -> double { return elem; }
  );
}

AlphabetEncoding<Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T> getEncoding(Rcpp::StringVector& input) {
  return getAlphabetEncoding<Rcpp::StringVector, Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T>(
      input,
      [](const Rcpp::String::StringProxy& elem) -> std::string { return Rcpp::as<std::string>(elem); }
  );
}