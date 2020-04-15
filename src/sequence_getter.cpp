#include "sequence_getter.h"

SequenceGetter_t<TidysqEncodedSequence> getTidysqRowGetter(std::vector<TidysqEncodedSequence>& encodedSequences) {
  return [&encodedSequences](int rowNum) -> TidysqEncodedSequence {
    return encodedSequences[rowNum];
  };
}