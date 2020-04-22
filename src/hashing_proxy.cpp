// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<vector>
#include "hash/polynomial_single_hasher.h"
#include "hash/complex_hasher.h"

//' @export
// [[Rcpp::export]]
int compute_polynomial_hash(int P,
                            int M,
                            Rcpp::IntegerVector items,
                            int begin,
                            int position) {
  PolynomialSingleHasher hasher(PolynomialSingleHasherConfig(P, M));
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

//' @export
// [[Rcpp::export]]
std::vector<int> compute_polynomial_multihash(Rcpp::IntegerVector P,
                                              Rcpp::IntegerVector M,
                                              Rcpp::IntegerVector items,
                                              int begin,
                                              int position) {
  std::vector<std::unique_ptr<SingleHasher>> singleHashers;
  for(int i = 0; i < P.size(); ++i) {
    singleHashers.emplace_back(
      new PolynomialSingleHasher(PolynomialSingleHasherConfig(P[i], M[i]))
    );
  }
  
  ComplexHasher hasher(std::move(singleHashers));
  for(const int& item: items) {
    hasher.append(item);
  }
  for(int i = 0; i < begin; ++i) {
    hasher.removeFirst(items[i]);
  }
  
  return position == -1 ?
    hasher.getHashes() :
    hasher.getHashes(position);
}
