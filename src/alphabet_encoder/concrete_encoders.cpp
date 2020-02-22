#include "alphabet_encoder.h"

std::unique_ptr<dictionary<int, int>> get_integer_alphabet_encoder() {
  return alphabet_encoder<int, int>(
      [](const int& elem) -> int { return elem; }
  )
}

