// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include "hash/primes.h"
#include "find_kmers_specific/find_kmers_integer_matrix.h"
#include "find_kmers_specific/find_kmers_string_matrix.h"
#include "find_kmers_specific/find_kmers_numeric_matrix.h"
#include "find_kmers_specific/find_kmers_string_list.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_contiguous_intervals_matrix(const Rcpp::IntegerVector &gaps) {
    auto intervals = gappedKMers::getContiguousIntervals(gaps);
    Rcpp::IntegerMatrix res(intervals.size(), 2);
    for (int i = 0; i < intervals.size(); ++i) {
        res(i, 0) = intervals[i].first + 1;
        res(i, 1) = intervals[i].second + 1;
    }
    return res;
}

inline std::vector<hashing::PolynomialSingleHasherConfig> getGappedKMerHasherConfigs(int hashDim) {
    std::vector<hashing::PolynomialSingleHasherConfig> res;
    for (int i = 0; i < hashDim; ++i) {
        res.emplace_back(hashing::hashPrimes[i].first, hashing::hashPrimes[i].second);
    }
    return std::move(res);
}

template<class sequences_t,
        class alphabet_t>
inline
Rcpp::List findKMers(sequences_t &sequences,
                     alphabet_t &alphabet,
                     std::vector<int> &gaps,
                     bool positionalKMers,
                     bool withKMerCounts,
                     const std::string &kmerDictionaryName,
                     int batchSize,
                     int hashDim,
                     bool verbose,
                     bool parallelMode) {
    auto hasherConfigs = std::move(getGappedKMerHasherConfigs(hashDim));
    return findKMersSpecific<decltype(hasherConfigs)>(
            sequences, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
            batchSize,
            verbose,
            parallelMode,
            hasherConfigs);
}

// [[Rcpp::export]]
Rcpp::List
find_gapped_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                         Rcpp::StringVector &alphabet,
                         std::vector<int> &gaps,
                         bool positionalKMers,
                         bool withKMerCounts,
                         const std::string &kmerDictionaryName,
                         int batchSize,
                         int hashDim,
                         bool verbose,
                         bool parallelMode) {
    return findKMers(sequenceMatrix, alphabet, gaps, positionalKMers, withKMerCounts,
                     kmerDictionaryName, batchSize, hashDim, verbose, parallelMode);
}

// [[Rcpp::export]]
Rcpp::List find_gapped_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                     Rcpp::IntegerVector &alphabet,
                                     std::vector<int> &gaps,
                                     bool positionalKMers,
                                     bool withKMerCounts,
                                     const std::string &kmerDictionaryName,
                                     int batchSize,
                                     int hashDim,
                                     bool verbose,
                                     bool parallelMode) {
    return findKMers(sequenceMatrix, alphabet, gaps, positionalKMers, withKMerCounts,
                     kmerDictionaryName, batchSize, hashDim, verbose, parallelMode);
}

// [[Rcpp::export]]
Rcpp::List find_gapped_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                     Rcpp::NumericVector &alphabet,
                                     std::vector<int> &gaps,
                                     bool positionalKMers,
                                     bool withKMerCounts,
                                     const std::string &kmerDictionaryName,
                                     int batchSize,
                                     int hashDim,
                                     bool verbose,
                                     bool parallelMode) {
    return findKMers(sequenceMatrix, alphabet, gaps, positionalKMers, withKMerCounts,
                     kmerDictionaryName, batchSize, hashDim, verbose, parallelMode);
}

// [[Rcpp::export]]
Rcpp::List find_gapped_kmers_list(Rcpp::List &sq,
                                  Rcpp::StringVector &alphabet,
                                  std::vector<int> &gaps,
                                  bool positionalKMers,
                                  bool withKMerCounts,
                                  const std::string &kmerDictionaryName,
                                  int batchSize,
                                  int hashDim,
                                  bool verbose,
                                  bool parallelMode) {
    return findKMers(sq, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, batchSize, hashDim,
                     verbose, parallelMode);
}
