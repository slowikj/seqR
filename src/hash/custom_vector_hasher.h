#ifndef SEQR_CUSTOM_VECTOR_HASHER_H
#define SEQR_CUSTOM_VECTOR_HASHER_H

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/functional/hash/hash.hpp>
#include <functional>
#include <algorithm>

namespace hashing::internal {

    template<class elem_t>
    inline std::size_t computeHash(const std::vector<elem_t> &v) {
        return boost::hash_range(v.begin(), v.end());
    }
}

namespace std {

    template<>
    struct hash<vector<int>> {
        size_t operator()(const vector<int> &c) const {
            return hashing::internal::computeHash(c);
        }
    };

    template<>
    struct hash<vector<long long>> {
        inline size_t operator()(const vector<long long> &c) const {
            return hashing::internal::computeHash(c);
        }
    };

    template<>
    struct hash<vector<uint64_t>> {
        inline size_t operator()(const vector<uint64_t> &c) const {
            return hashing::internal::computeHash(c);
        }
    };

}

#endif //SEQR_CUSTOM_VECTOR_HASHER_H
