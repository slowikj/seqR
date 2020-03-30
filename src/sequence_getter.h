#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <functional>
#include <Rcpp.h>
#include <memory>

template <class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

template <class input_matrix_t, class input_vector_t>
SequenceGetter_t<input_vector_t> getRcppMatrixRowGetter(input_matrix_t& sequenceMatrix) {
  return [&sequenceMatrix](int rowNum) {
    return std::move(sequenceMatrix(rowNum, Rcpp::_));
  };
}

#endif