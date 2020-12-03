// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <vector>
#include "hash/primes.h"
#include "find_kmers_specific/find_kmers_integer_matrix.h"
#include "find_kmers_specific/find_kmers_string_matrix.h"
#include "find_kmers_specific/find_kmers_numeric_matrix.h"
#include "find_kmers_specific/find_kmers_string_list.h"
#include "hash/polynomial_single_hasher.h"

inline hashing::ComplexHasher createKMerComplexHasher(int hashDim) {
    std::vector<std::unique_ptr<hashing::SingleHasher>> singleHashers;
    for (int i = 0; i < hashDim; ++i) {
        singleHashers.emplace_back(
                new hashing::PolynomialSingleHasher(
                        hashing::PolynomialSingleHasherConfig(hashing::hashPrimes[i].first,
                                                              hashing::hashPrimes[i].second)));
    }
    hashing::ComplexHasher complexHasher(std::move(singleHashers));
    return complexHasher;
}

template<class sequences_t,
        class alphabet_t>
Rcpp::List countContiguousKMers(
        sequences_t &sequences,
        alphabet_t &alphabet,
        int k,
        bool positionalKMers,
        bool withKMerCounts,
        const std::string &kmerDictionaryName,
        int batchSize,
        int hashDim,
        bool verbose,
        bool parallelMode) {
    std::function<hashing::ComplexHasher()> algorithmParams = [hashDim]() -> hashing::ComplexHasher {
        return createKMerComplexHasher(hashDim);
    };
    std::vector<int> gaps(k - 1);
    return findKMersSpecific<decltype(algorithmParams)>(
            sequences, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
            batchSize,
            verbose,
            parallelMode,
            algorithmParams);
}

// [[Rcpp::export]]
Rcpp::List count_contiguous_kmers_string(
        Rcpp::StringMatrix &sequenceMatrix,
        Rcpp::StringVector &alphabet,
        int k,
        bool positionalKMers,
        bool withKMerCounts,
        const std::string &kmerDictionaryName,
        int batchSize,
        int hashDim,
        bool verbose,
        bool parallelMode) {
    return countContiguousKMers(
            sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts,
            kmerDictionaryName,
            batchSize,
            hashDim,
            verbose,
            parallelMode
    );
}

// [[Rcpp::export]]
Rcpp::List count_contiguous_kmers_integer(
        Rcpp::IntegerMatrix &sequenceMatrix,
        Rcpp::IntegerVector &alphabet,
        int k,
        bool positionalKMers,
        bool withKMerCounts,
        const std::string &kmerDictionaryName,
        int batchSize,
        int hashDim,
        bool verbose,
        bool parallelMode) {
    return countContiguousKMers(
            sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts,
            kmerDictionaryName,
            batchSize,
            hashDim,
            verbose,
            parallelMode
    );
}

// [[Rcpp::export]]
Rcpp::List count_contiguous_kmers_numeric(
        Rcpp::NumericMatrix &sequenceMatrix,
        Rcpp::NumericVector &alphabet,
        int k,
        bool positionalKMers,
        bool withKMerCounts,
        const std::string &kmerDictionaryName,
        int batchSize,
        int hashDim,
        bool verbose,
        bool parallelMode) {
    return countContiguousKMers(
            sequenceMatrix, alphabet, k, positionalKMers, withKMerCounts,
            kmerDictionaryName,
            batchSize,
            hashDim,
            verbose,
            parallelMode
    );
}

// [[Rcpp::export]]
Rcpp::List count_contiguous_kmers_list(
        Rcpp::List &sq,
        Rcpp::StringVector &alphabet,
        int k,
        bool positionalKMers,
        bool withKMerCounts,
        const std::string &kmerDictionaryName,
        int batchSize,
        int hashDim,
        bool verbose,
        bool parallelMode) {
    return countContiguousKMers(
            sq, alphabet, k, positionalKMers, withKMerCounts, kmerDictionaryName, batchSize,
            hashDim,
            verbose,
            parallelMode);
}
