#pragma once

// [[Rcpp::depends(BH)]]
#include <boost/functional/hash/hash.hpp>
#include <functional>
#include <algorithm>

namespace hashing::internal
{

    template <class elem_t>
    inline std::size_t computeHash(const std::vector<elem_t> &v)
    {
        return boost::hash_range(v.begin(), v.end());
    }
}

namespace std
{

    template <>
    struct hash<vector<uint32_t>>
    {
        size_t operator()(const vector<uint32_t> &c) const
        {
            return hashing::internal::computeHash(c);
        }
    };

    template <>
    struct hash<vector<uint64_t>>
    {
        inline size_t operator()(const vector<uint64_t> &c) const
        {
            return hashing::internal::computeHash(c);
        }
    };

}
