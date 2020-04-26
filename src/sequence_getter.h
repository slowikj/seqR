#ifndef SEQUENCE_GETTER_H
#define SEQUENCE_GETTER_H

#include <Rcpp.h>
#include <functional>
#include <memory>
#include <vector>
#include "safe_sequences_wrapper.h"

template<class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

inline SequenceGetter_t<SafeTidysqSequencesWrapper::Row>
getTidysqRowGetter(SafeTidysqSequencesWrapper& safeWrapper) {
    return [&safeWrapper](int rowNum) -> SafeTidysqSequencesWrapper::Row {
        return std::move(safeWrapper.row(rowNum));
    };
}

template<class elem_t>
inline SequenceGetter_t<typename SafeSequencesWrapper<elem_t>::Row>
getSafeMatrixRowGetter(SafeMatrixSequenceWrapper<elem_t>& sequenceWrapper) {
    return [&sequenceWrapper](int rowNum) -> typename SafeSequencesWrapper<elem_t>::Row {
        return std::move(sequenceWrapper.row(rowNum));
    };
}

#endif
