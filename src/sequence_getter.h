#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <Rcpp.h>
#include <functional>
#include <memory>
#include <vector>

template <class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

SequenceGetter_t<Rcpp::Fast<Rcpp::RawVector>> getTidysqRowGetter(std::vector<Rcpp::RawVector>& encodedSequences);

template <class input_matrix_t, class input_vector_t>
SequenceGetter_t<input_vector_t> getRcppMatrixRowGetter(input_matrix_t& sequenceMatrix) {
  return [&sequenceMatrix](int rowNum) -> input_vector_t {
    return sequenceMatrix(rowNum, Rcpp::_);
  };
}

#endif
