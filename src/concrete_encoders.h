#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

#include <Rcpp.h>

#include "alphabet_encoder.h"
#include <memory>
#include <string>

typedef short ENCODED_ELEM_T;

AlphabetEncoding<int, int, ENCODED_ELEM_T> getEncoding(Rcpp::IntegerVector& input);
  
AlphabetEncoding<double, double, ENCODED_ELEM_T> getEncoding(Rcpp::NumericVector& input);

AlphabetEncoding<Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T> getEncoding(Rcpp::StringVector& input);

#endif //CONCRETE_ENCODERS_H
