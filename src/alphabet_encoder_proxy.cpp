// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<iostream>

#include "alphabet_encoder/concrete_encoders.h"
#include <utility>

template<class rcpp_input_t>
Rcpp::IntegerVector prepareOutputVector(rcpp_input_t& input) {
  auto alphabetEncoding = getEncoding(input);
  auto res = Rcpp::IntegerVector(input.size());
  for(int i = 0; i < input.size(); ++i) {
    res[i] = alphabetEncoding.encode(input(i));
  }
  res.names() = input;
  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_integer_alphabet(Rcpp::IntegerVector& input) {
  return prepareOutputVector<Rcpp::IntegerVector>(input);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_numeric_alphabet(Rcpp::NumericVector& input) {
  return prepareOutputVector<Rcpp::NumericVector>(input);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_string_alphabet(Rcpp::StringVector& input) {
  return prepareOutputVector<Rcpp::StringVector>(input);
}
