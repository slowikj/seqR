// [[Rcpp::plugins("cpp17")]]
#include<Rcpp.h>
#include<vector>
#include "hash/polynomial_single_hasher.h"
#include "hash/complex_hasher.h"

//' @export
// [[Rcpp::export]]
int compute_polynomial_hash(int P,
                            int M,
                            Rcpp::IntegerVector items,
                            int begin) {
    hashing::PolynomialSingleHasher hasher(hashing::PolynomialSingleHasherConfig(P, M));
    for (const int &item: items) {
        hasher.append(item);
    }
    for (int i = 0; i < begin; ++i) {
        hasher.removeFirst(items[i]);
    }
    return hasher.getHash();
}

//' @export
// [[Rcpp::export]]
std::vector<int> compute_polynomial_multihash(Rcpp::IntegerVector P,
                                              Rcpp::IntegerVector M,
                                              Rcpp::IntegerVector items,
                                              int begin,
                                              int position) {
    std::vector<std::unique_ptr<hashing::SingleHasher>> singleHashers;
    for (int i = 0; i < P.size(); ++i) {
        singleHashers.emplace_back(
                new hashing::PolynomialSingleHasher(hashing::PolynomialSingleHasherConfig(P[i], M[i]))
        );
    }

    hashing::ComplexHasher hasher(std::move(singleHashers));
    for (const int &item: items) {
        hasher.append(item);
    }
    for (int i = 0; i < begin; ++i) {
        hasher.removeFirst(items[i]);
    }

    return position == -1 ?
           hasher.getHashes() :
           hasher.getHashes(position - 1);
}
