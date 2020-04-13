// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>
#include "alphabet_encoder.h"
#include "app_conf.h"
#include "hash/custom_hashers.h"
#include <utility>
#include <functional>

template <class rcpp_input_t, class input_elem_t, class encoded_elem_t, class hasher_t>
Rcpp::IntegerVector prepareOutputVector(
    rcpp_input_t& input) {
  auto alphabetEncoding = getAlphabetEncoding<rcpp_input_t, input_elem_t, encoded_elem_t, hasher_t>(
    input
  );
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
  return prepareOutputVector<Rcpp::IntegerVector,
                             int,
                             ENCODED_ELEM_T,
                             std::hash<int>>(input);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_numeric_alphabet(Rcpp::NumericVector& input) {
  return prepareOutputVector<Rcpp::NumericVector,
                             double,
                             ENCODED_ELEM_T,
                             std::hash<double>>(input);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_string_alphabet(Rcpp::StringVector& input) {
  return prepareOutputVector<Rcpp::StringVector,
                             Rcpp::StringVector::stored_type,
                             ENCODED_ELEM_T,
                             string_proxy_hasher>(input);
}
