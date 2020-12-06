#ifndef SEQR_TYPES_H
#define SEQR_TYPES_H

#include "custom_vector_hasher.h"
#include <vector>

namespace hashing {
    using single_hash_t = long long;

    using multidim_hash_t = std::vector<single_hash_t>;

    using multidim_hasher_t = LLVectorHasher;
}

#endif //SEQR_TYPES_H
