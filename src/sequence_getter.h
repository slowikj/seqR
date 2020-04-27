#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <Rcpp.h>
#include <functional>
#include <memory>
#include <vector>
#include "safe_sequences_wrapper.h"

template<class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

template<class elem_t>
inline SequenceGetter_t<typename SafeSequencesWrapper<elem_t>::Row>
getSafeMatrixRowGetter(SafeMatrixSequenceWrapper<elem_t> &sequenceWrapper) {
    return [&sequenceWrapper](int rowNum) -> typename SafeSequencesWrapper<elem_t>::Row {
        return std::move(sequenceWrapper.row(rowNum));
    };
}

template<class rcpp_matrix_t, class elem_t>
inline SequenceGetter_t<typename RcppParallel::RMatrix<elem_t>::Row>
getRMatrixRowGetter(rcpp_matrix_t &rcppMatrix) {
    RcppParallel::RMatrix<elem_t> wrappedMatrix(rcppMatrix);
    return [wrappedMatrix](int rowNum) -> typename RcppParallel::RMatrix<elem_t>::Row {
        return wrappedMatrix.row(rowNum);
    };
}

inline SequenceGetter_t<SafeTidysqSequencesWrapper::Row>
getTidysqRowGetter(SafeTidysqSequencesWrapper &safeWrapper) {
    return [&safeWrapper](int rowNum) -> SafeTidysqSequencesWrapper::Row {
        return std::move(safeWrapper.row(rowNum));
    };
}

inline SequenceGetter_t<RcppParallel::RVector<unsigned char>>
getTidysqRVectorGetter(Rcpp::List &unpackedSequences) {
    return [&unpackedSequences](int rowNum) -> RcppParallel::RVector<unsigned char> {
        Rcpp::RawVector sequence = unpackedSequences[rowNum];
        return RcppParallel::RVector<unsigned char>(sequence);
    };
}

#endif
