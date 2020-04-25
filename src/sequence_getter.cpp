// [[Rcpp::plugins("cpp17")]]
#include "sequence_getter.h"

SequenceGetter_t<SafeTidysqSequencesWrapper::Row> getTidysqRowGetter(SafeTidysqSequencesWrapper& safeWrapper) {
    return [&safeWrapper](int rowNum) -> SafeTidysqSequencesWrapper::Row {
        return std::move(safeWrapper.row(rowNum));
    };
}
