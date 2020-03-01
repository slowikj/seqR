// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<iostream>

#include "alphabet_encoder/concrete_encoders.h"
#include <utility>

template<class rcpp_input_t, class rcpp_item_t>
rcpp_input_t prepareOutputVector(const rcpp_input_t& input) {
  auto alphabetEncoding = getEncoding(input);
  auto res = rcpp_input_t(alphabetEncoding.alphabetSize());
  for(int i = 0; i < input.size(); ++i) {
    res[i] = alphabetEncoding.encode(input(i));
  }
  res.names() = input;
  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_integer_alphabet(const Rcpp::IntegerVector& input) {
  return prepareOutputVector<Rcpp::IntegerVector, int>(input);
}
