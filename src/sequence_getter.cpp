#include "sequence_getter.h"

SequenceGetter_t<Rcpp::RawVector> getTidysqRowGetter(std::vector<Rcpp::RawVector>& encodedSequences) {
  return [&encodedSequences](int rowNum) -> Rcpp::RawVector {
    return std::move(encodedSequences[rowNum]);
  };
}
