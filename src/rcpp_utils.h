#ifndef RCPP_UTILS_H
#define RCPP_UTILS_H

#include <Rcpp.h>
#include <memory>
#include "kmer_counting_common.h"

template <class input_matrix_t, class input_vector_t>
RowGetter_t<input_vector_t> getRcppMatrixRowGetter(input_matrix_t& sequenceMatrix) {
  return [&sequenceMatrix](int rowNum) {
    return std::move(sequenceMatrix(rowNum, Rcpp::_));
  };
}

#endif