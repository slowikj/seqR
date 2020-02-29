#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>

#include "alphabet_encoder.h"
#include "../dictionary.h"
#include <memory>

// INTEGER ENCODING
AlphabetEncoder<int, int> getIntegerAlphabetEncoder() {
  return AlphabetEncoder<int, int>([](const int& elem) -> int { return elem; });
}

Dictionary<int, short> getEncoding(const Rcpp::IntegerVector& input) {
  auto encoder = getIntegerAlphabetEncoder();
  return encoder.getEncoding<short, Rcpp::IntegerVector>(input);
}

#endif //CONCRETE_ENCODERS_H