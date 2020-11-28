// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include "default_alphabet_encoder.h"
#include "hash/custom_hashers.h"
#include "dictionary/unordered_map_wrapper.h"
#include <utility>
#include <functional>

template<class rcpp_input_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
Rcpp::IntegerVector prepareOutputVector(
        rcpp_input_t &input) {
    auto alphabetEncoding = getDefaultAlphabetEncoder<rcpp_input_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
            input
    );
    auto res = Rcpp::IntegerVector(input.size());
    for (int i = 0; i < input.size(); ++i) {
        res[i] = alphabetEncoding.encode(input(i));
    }
    res.names() = input;
    return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_integer_alphabet(Rcpp::IntegerVector &input) {
    return prepareOutputVector<Rcpp::IntegerVector,
            int,
            short,
            UnorderedMapWrapper>(input);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_numeric_alphabet(Rcpp::NumericVector &input) {
    return prepareOutputVector<Rcpp::NumericVector,
            double,
            short,
            UnorderedMapWrapper>(input);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector encode_string_alphabet(Rcpp::StringVector &input) {
    return prepareOutputVector<Rcpp::StringVector,
            Rcpp::StringVector::stored_type,
            short,
            UnorderedMapWrapper>(input);
}
