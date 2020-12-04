// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include "hash/primes.h"
#include "count_kmers_specific/count_kmers_integer_matrix.h"
#include "count_kmers_specific/count_kmers_string_matrix.h"
#include "count_kmers_specific/count_kmers_numeric_matrix.h"
#include "count_kmers_specific/count_kmers_string_list.h"

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
Rcpp::List countGappedKMers(
        sequences_t &sequences,
        alphabet_t &alphabet,
        Rcpp::Environment &rcppParams) {
    auto userParams = std::move(UserParams::createForGapped(rcppParams));
    auto hasherConfigs = std::move(getGappedKMerHasherConfigs(userParams.hashDim));
    return countKMersSpecific<decltype(hasherConfigs)>(
            sequences, alphabet, userParams, hasherConfigs);
}

// [[Rcpp::export(".count_gapped_kmers_string")]]
Rcpp::List count_gapped_kmers_string(
        Rcpp::StringMatrix &sequenceMatrix,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(
            sequenceMatrix, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_integer")]]
Rcpp::List count_gapped_kmers_integer(
        Rcpp::IntegerMatrix &sequenceMatrix,
        Rcpp::IntegerVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(
            sequenceMatrix, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_numeric")]]
Rcpp::List count_gapped_kmers_numeric(
        Rcpp::NumericMatrix &sequenceMatrix,
        Rcpp::NumericVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(
            sequenceMatrix, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_list")]]
Rcpp::List count_gapped_kmers_list(
        Rcpp::List &sq,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(
            sq, alphabet, rcppParams);
}
