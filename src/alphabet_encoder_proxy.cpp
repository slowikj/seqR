// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>
#include "input_to_internal_item_converter.h"
#include "alphabet_encoder.h"
#include "app_conf.h"
#include <utility>

template <class rcpp_input_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
Rcpp::IntegerVector prepareOutputVector(
    rcpp_input_t& input,
    InputToInternalItemConverter_t<input_elem_t, internal_elem_t> inputToInternalItemConverter) {
  auto alphabetEncoding = getAlphabetEncoding<rcpp_input_t, input_elem_t, internal_elem_t, encoded_elem_t>(
    input,
    inputToInternalItemConverter
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
  return prepareOutputVector<Rcpp::IntegerVector, int, int, ENCODED_ELEM_T>(
      input,
      getIntToIntConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_numeric_alphabet(Rcpp::NumericVector& input) {
  return prepareOutputVector<Rcpp::NumericVector, double, double, ENCODED_ELEM_T>(
      input,
      getDoubleToDoubleConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_string_alphabet(Rcpp::StringVector& input) {
  return prepareOutputVector<Rcpp::StringVector, Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T>(
      input,
      getRcppStringToStringConverter());
}
