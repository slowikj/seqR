#ifndef CUSTOM_HASHERS_H
#define CUSTOM_HASHERS_H

#include <Rcpp.h>
#include <functional>
#include <algorithm>

template<>
struct std::hash<Rcpp::StringVector::stored_type> {
    inline std::size_t operator()(const Rcpp::StringVector::stored_type &v) const {
        return rcppStringHasher(v);
    }

private:
    std::hash<Rcpp::String> rcppStringHasher;
};



#endif
