// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>

#include "hash/polynomial_single_hasher.h"

//' @export
// [[Rcpp::export]]
int compute_polynomial_hash(int P,
                            int M,
                            Rcpp::IntegerVector items,
                            int begin,
                            int position) {
  PolynomialSingleHasher hasher(P, M);
  for(const int& item: items) {
    hasher.append(item);
  }
  for(int i = 0; i < begin; ++i) {
    hasher.removeFirst(items[i]);
  }
  return position == -1 ?
    hasher.getHash() :
    hasher.getHash(position);
}
