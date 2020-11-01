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
Rcpp::List findKMers(sequences_t &sequences,
                     alphabet_t &alphabet,
                     int k,
                     bool positionalKMers,
                     bool withKMerCounts,
                     const std::string &kmerDictionaryName,
                     int batchSize) {
    std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
    std::vector<int> gaps(k - 1);
    return findKMersSpecific<decltype(algorithmParams)>(
            sequences, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
            batchSize,
            algorithmParams);
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                             std::vector<std::string> &alphabet,
                             int k,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize) {
    return findKMers(sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts,
                     kmerDictionaryName,
                     batchSize
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                              std::vector<int> &alphabet,
                              int k,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName,
                              int batchSize) {
    return findKMers(sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts,
                     kmerDictionaryName,
                     batchSize
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                              std::vector<double> &alphabet,
                              int k,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName,
                              int batchSize) {
    return findKMers(sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts,
                     kmerDictionaryName,
                     batchSize
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_kmers_list(Rcpp::List &sq,
                           Rcpp::StringVector &alphabet,
                           int k,
                           bool positionalKMers,
                           bool withKMerCounts,
                           const std::string &kmerDictionaryName,
                           int batchSize) {
    return findKMers(sq, alphabet, k, positionalKMers, withKMerCounts, kmerDictionaryName, batchSize);
}
