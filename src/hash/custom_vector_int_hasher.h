#ifndef SEQR_CUSTOM_VECTOR_INT_HASHER_H
#define SEQR_CUSTOM_VECTOR_INT_HASHER_H

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/functional/hash/hash.hpp>
#include <functional>
#include <algorithm>

namespace hashing {

    struct IntVectorHasher {
        inline std::size_t operator()(const std::vector<int> &c) const {
            return boost::hash_range(c.begin(), c.end());
        }
    };
}

#endif //SEQR_CUSTOM_VECTOR_INT_HASHER_H
