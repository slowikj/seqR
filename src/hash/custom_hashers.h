#ifndef CUSTOM_HASHERS_H
#define CUSTOM_HASHERS_H

#include <Rcpp.h>
#include <functional>

// we assume that a given StringVector stores only single characters
struct string_proxy_hasher {
  
  std::size_t operator()(const Rcpp::StringVector::stored_type& v) const {
    // TODO: improve the function
    return std::hash<std::string>()(std::move(Rcpp::as<std::string>(v)));
  }
};

struct vector_int_hasher {
  
  const static int P = 47;
  
  const static int M = 1e9 + 7;
  
  std::size_t operator()(const std::vector<int> &c) const {
    // TODO: try to use boost::hash_range(c.begin(), c.end());
    return std::accumulate(std::begin(c), std::end(c),
                           0,
                           [](size_t prev, int elem) -> int {
                             return static_cast<int>((static_cast<long long>(prev) * P + elem) % M);
                           });
  }
};

#endif
