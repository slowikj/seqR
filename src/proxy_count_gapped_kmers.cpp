// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include "hash/globals.h"
#include "kmer_task_param_dispatcher.h"

inline std::vector<hashing::PolynomialSingleHasherConfig> getGappedKMerHasherConfigs(int hashDim) {
    std::vector<hashing::PolynomialSingleHasherConfig> res;
    for (int i = 0; i < hashDim; ++i) {
        res.emplace_back(hashing::config::hashPrimes[i].first, hashing::config::hashPrimes[i].second);
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
    return countKMers<sequences_t, alphabet_t, decltype(hasherConfigs)>(
            sequences, alphabet, userParams, hasherConfigs);
}

// [[Rcpp::export(".count_gapped_kmers_string")]]
Rcpp::List count_gapped_kmers_string(
        Rcpp::StringMatrix &sequenceMatrix,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(sequenceMatrix, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_integer")]]
Rcpp::List count_gapped_kmers_integer(
        Rcpp::IntegerMatrix &sequenceMatrix,
        Rcpp::IntegerVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(sequenceMatrix, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_numeric")]]
Rcpp::List count_gapped_kmers_numeric(
        Rcpp::NumericMatrix &sequenceMatrix,
        Rcpp::NumericVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(sequenceMatrix, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_string_vector")]]
Rcpp::List count_gapped_kmers_string_vector(
        Rcpp::StringVector &sq,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(sq, alphabet, rcppParams);
}

// [[Rcpp::export(".count_gapped_kmers_tidysq")]]
Rcpp::List count_gapped_kmers_tidysq(
        Rcpp::List &sq,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countGappedKMers(sq, alphabet, rcppParams);
}
