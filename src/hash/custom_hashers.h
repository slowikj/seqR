#ifndef CUSTOM_HASHERS_H
#define CUSTOM_HASHERS_H

#include <Rcpp.h>

// we assume that a given StringVector stores only single characters
struct string_proxy_hasher {
  
  const int P = 107;
  
  const int M = 1e9 + 9; 
  
  std::size_t operator()(const Rcpp::StringVector::stored_type& v) const {
    // TODO: improve the function
    std::string str = Rcpp::as<std::string>(v);
    long long res = 0;
    for(const auto& elem: str) {
      res = ((res * P) + elem) % M;
    }
    return static_cast<std::size_t>(res);
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
