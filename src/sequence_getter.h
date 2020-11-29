#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include<RcppParallel.h>
#include <functional>
#include <memory>
#include <vector>
#include "safe_sequences_wrapper.h"

template<class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

template<class rcpp_matrix_t, class elem_t>
inline SequenceGetter_t<typename RcppParallel::RMatrix<elem_t>::Row>
getRMatrixRowGetter(rcpp_matrix_t &rcppMatrix, int rowOffset = 0) {
    RcppParallel::RMatrix<elem_t> wrappedMatrix(rcppMatrix);
    return [wrappedMatrix, rowOffset](int rowNum) -> typename RcppParallel::RMatrix<elem_t>::Row {
        return wrappedMatrix.row(rowNum + rowOffset);
    };
}



#endif
