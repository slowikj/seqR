#ifndef CONCRETE_ENCODERS_H
#define CONCRETE_ENCODERS_H

#include "alphabet_encoder.h"
#include <memory>

std::unique_ptr<alphabet_encoder<int, int>> get_integer_alphabet_encoder();

#endif //CONCRETE_ENCODERS_H