// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <vector>
#include <functional>
#include "kmer_task_solver.h"

inline ComplexHasher createKMerComplexHasher() {
    std::vector<std::unique_ptr<SingleHasher>> singleHashers;
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(101, 1e9 + 33)));
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(97, 1e9 + 7)));
    ComplexHasher complexHasher(std::move(singleHashers));
    return complexHasher;
}

template<class sequences_t,
        class alphabet_t>
Rcpp::IntegerMatrix findKMers(sequences_t &sequences,
                              alphabet_t &alphabet,
                              int k,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName) {
    std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
    std::vector<int> gaps(k - 1);
    return findKMers<decltype(algorithmParams)>(
            sequences, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, algorithmParams);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                                      Rcpp::StringVector &alphabet,
                                      int k,
                                      bool positionalKMers,
                                      bool withKMerCounts,
                                      const std::string &kmerDictionaryName) {
    return findKMers<Rcpp::StringMatrix, Rcpp::StringVector>(
            sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts, kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                       Rcpp::IntegerVector &alphabet,
                                       int k,
                                       bool positionalKMers,
                                       bool withKMerCounts,
                                       const std::string &kmerDictionaryName) {
    return findKMers<Rcpp::IntegerMatrix, Rcpp::IntegerVector>(
            sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts, kmerDictionaryName
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                       Rcpp::NumericVector &alphabet,
                                       int k,
                                       bool positionalKMers,
                                       bool withKMerCounts,
                                       const std::string &kmerDictionaryName) {
    return findKMers<Rcpp::NumericMatrix, Rcpp::NumericVector>(
            sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts, kmerDictionaryName
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_tidysq(Rcpp::List &sq,
                                      Rcpp::StringVector &alphabet,
                                      int k,
                                      bool positionalKMers,
                                      bool withKMerCounts,
                                      const std::string &kmerDictionaryName) {
    return findKMers<Rcpp::List, Rcpp::StringVector>(
            sq, alphabet, k, positionalKMers, withKMerCounts, kmerDictionaryName
    );
}
