#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>

#include "alphabet_encoder.h"
#include "../dictionary.h"
#include <memory>

typedef int ENCODED_T;

// INTEGER ENCODING
AlphabetEncoder<int, int> getIntegerAlphabetEncoder() {
  return AlphabetEncoder<int, int>([](const int& elem) -> int { return elem; });
}

Dictionary<int, ENCODED_T> getEncoding(const Rcpp::IntegerVector& input) {
  auto encoder = getIntegerAlphabetEncoder();
  return encoder.getEncoding<ENCODED_T, Rcpp::IntegerVector>(input);
}

#endif //CONCRETE_ENCODERS_H