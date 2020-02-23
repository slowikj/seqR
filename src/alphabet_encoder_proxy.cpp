// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>

#include "alphabet_encoder/concrete_encoders.h"

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector encode_alphabet(const Rcpp::NumericVector& input) {
  auto integer_encoder = get_integer_alphabet_encoder();
  auto encoded = integer_encoder.get_encoding<short, Rcpp::NumericVector>(input);
  for(const auto& elem: encoded) {
    Rcpp::Rcout << elem.first << " " << elem.second << std::endl;
  }
  
  return Rcpp::NumericVector::create(1,2,3);
}