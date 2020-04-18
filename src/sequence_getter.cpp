#include "sequence_getter.h"

SequenceGetter_t<TidysqEncodedSequenceProxy> getTidysqRowGetter(
    std::vector<TidysqEncodedSequenceProxy>& encodedSequences) {
  return [&encodedSequences](int rowNum) -> TidysqEncodedSequenceProxy {
    return encodedSequences[rowNum];
  };
}
