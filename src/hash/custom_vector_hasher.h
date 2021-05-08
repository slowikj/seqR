#pragma once

#include <functional>
#include <algorithm>

namespace hashing::internal
{

    template <class elem_t>
    inline std::size_t computeHash(const std::vector<elem_t> &v)
    {
        // implementation: boost::hash_range(v.begin(), v.end());
        std::size_t res = 0;
        for(const auto& elem: v)
        {
            res ^= elem + 0x9e3779b9 + (res << 6) + (res >> 2);
        }
        return res;
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
