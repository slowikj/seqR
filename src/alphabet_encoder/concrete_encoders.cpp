#include "alphabet_encoder.h"

AlphabetEncoder<int, int> get_integer_alphabet_encoder() {
  return AlphabetEncoder<int, int>([](const int& elem) -> int { return elem; });
}

