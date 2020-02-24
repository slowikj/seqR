// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<iostream>

#include "alphabet_encoder/concrete_encoders.h"
#include "alphabet_encoder/alphabet_encoder.h"
#include <utility>

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_alphabet(Rcpp::IntegerVector input) {
  AlphabetEncoder<int,int> integer_encoder = getIntegerAlphabetEncoder();
  Dictionary<int, short> encoded = integer_encoder.getEncoding<short, Rcpp::IntegerVector>(input);
  
  auto res = Rcpp::IntegerVector(encoded.size());
  auto keys = encoded.getKeys();
  
  for(int i = 0; i < keys.size(); ++i) {
    res[i] = encoded[keys[i]];
  }
  res.names() = keys;
  
  return res;
}