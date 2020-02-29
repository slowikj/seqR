// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<iostream>

#include "alphabet_encoder/concrete_encoders.h"
#include <utility>

template<class rcpp_item_t, class rcpp_input_t>
rcpp_input_t prepareOutputVector(const rcpp_input_t& input) {
  Dictionary<rcpp_item_t, ENCODED_T> encoded = getEncoding(input);
  auto res = Rcpp::IntegerVector(encoded.size());
  auto keys = encoded.getKeys();
  for(int i = 0; i < keys.size(); ++i) {
    res[i] = encoded[keys[i]];
  }
  res.names() = keys;
  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_integer_alphabet(Rcpp::IntegerVector input) {
  return prepareOutputVector<int, Rcpp::IntegerVector>(input);
}