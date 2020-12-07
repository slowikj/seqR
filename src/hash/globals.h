#ifndef SEQR_GLOBALS_H
#define SEQR_GLOBALS_H

#include "custom_vector_hasher.h"
#include <vector>

namespace hashing::config {

    using single_hash_t = long long;

    using multidim_hash_t = std::vector<single_hash_t>;

    using multidim_hasher_t = LLVectorHasher;

    const static std::vector<std::pair<int, int>> hashPrimes = {
            {101, 1e9 + 7},
            {97,  1e9 + 33},
            {89,  1e9 + 9},
            {71,  1192901131},
            {67,  1190705431},
            {61,  1190494969},
            {59,  1086218491},
            {53,  1065433423}
    };
}

#endif //SEQR_GLOBALS_H
