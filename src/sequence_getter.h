#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <Rcpp.h>
#include <functional>
#include <memory>
#include <vector>
#include "tidysq_encoded_sequence.h"

template <class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

SequenceGetter_t<TidysqEncodedSequence> getTidysqRowGetter(std::vector<TidysqEncodedSequence>& encodedSequences);

template <class input_matrix_t, class input_vector_t>
SequenceGetter_t<input_vector_t> getRcppMatrixRowGetter(input_matrix_t& sequenceMatrix) {
  return [&sequenceMatrix](int rowNum) -> input_vector_t {
    return std::move(sequenceMatrix(rowNum, Rcpp::_));
  };
}

#endif