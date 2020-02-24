// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<iostream>

#include "alphabet_encoder/concrete_encoders.h"
#include "alphabet_encoder/alphabet_encoder.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_alphabet(Rcpp::IntegerVector input) {
  AlphabetEncoder<int,int> integer_encoder = getIntegerAlphabetEncoder();
  Dictionary<int, short> encoded = integer_encoder.getEncoding<short, Rcpp::IntegerVector>(input);
  for(const auto& elem: encoded) {
    Rcpp::Rcout << elem.first << " " << elem.second << std::endl;
  }
  
  return Rcpp::IntegerVector::create(1,2,3);
}