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

template<class elem_t>
inline SequenceGetter_t<typename SafeSequencesWrapper<elem_t>::Row>
getSafeMatrixRowGetter(SafeMatrixSequenceWrapper<elem_t> &sequenceWrapper, int offset = 0) {
    return [&sequenceWrapper, offset](int rowNum) -> typename SafeSequencesWrapper<elem_t>::Row {
        return std::move(sequenceWrapper.row(rowNum + offset));
    };
}

template<class rcpp_matrix_t, class elem_t>
inline SequenceGetter_t<typename RcppParallel::RMatrix<elem_t>::Row>
getRMatrixRowGetter(rcpp_matrix_t &rcppMatrix, int offset = 0) {
    RcppParallel::RMatrix<elem_t> wrappedMatrix(rcppMatrix);
    return [wrappedMatrix, offset](int rowNum) -> typename RcppParallel::RMatrix<elem_t>::Row {
        return wrappedMatrix.row(rowNum + offset);
    };
}

inline SequenceGetter_t<SafeTidysqSequencesWrapper::Row>
getTidysqRowGetter(SafeTidysqSequencesWrapper &safeWrapper, int offset = 0) {
    return [&safeWrapper, offset](int rowNum) -> SafeTidysqSequencesWrapper::Row {
        return std::move(safeWrapper.row(rowNum + offset));
    };
}

inline SequenceGetter_t<RcppParallel::RVector<unsigned char>>
getTidysqRVectorGetter(std::vector<RcppParallel::RVector<unsigned char>> &unpackedSafeSequences, int offset = 0) {
    return [&unpackedSafeSequences, offset](int rowNum) -> RcppParallel::RVector<unsigned char> {
        return unpackedSafeSequences[rowNum + offset];
    };
}

#endif
