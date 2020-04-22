#include "sequence_getter.h"

SequenceGetter_t<Rcpp::Fast<Rcpp::RawVector>> getTidysqRowGetter(std::vector<Rcpp::RawVector>& encodedSequences) {
  return [&encodedSequences](int rowNum) -> Rcpp::Fast<Rcpp::RawVector> {
    return Rcpp::Fast<Rcpp::RawVector>(encodedSequences[rowNum]);
  };
}
