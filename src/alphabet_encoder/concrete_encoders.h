#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

#include "alphabet_encoder.h"
#include <memory>

AlphabetEncoder<int, int> getIntegerAlphabetEncoder() {
  return AlphabetEncoder<int, int>([](const int& elem) -> int { return elem; });
}

#endif //CONCRETE_ENCODERS_H