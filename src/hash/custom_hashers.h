#ifndef CUSTOM_HASHERS_H
#define CUSTOM_HASHERS_H

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/functional/hash/hash.hpp>
#include <functional>
#include <algorithm>

struct string_proxy_hasher {
    inline std::size_t operator()(const Rcpp::StringVector::stored_type &v) const {
        return rcppStringHasher(v);
    }

private:
    std::hash<Rcpp::String> rcppStringHasher;
};

template<>
struct std::hash<std::vector<int>> {

    const static int P = 47;

    const static int M = 1e9 + 7;

    inline std::size_t operator()(const std::vector<int> &c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

#endif
