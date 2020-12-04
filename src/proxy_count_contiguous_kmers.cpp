// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <vector>
#include "hash/primes.h"
#include "count_kmers_specific/count_kmers_integer_matrix.h"
#include "count_kmers_specific/count_kmers_string_matrix.h"
#include "count_kmers_specific/count_kmers_numeric_matrix.h"
#include "count_kmers_specific/count_kmers_string_list.h"
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
        Rcpp::Environment &rcppParams) {
    auto userParams = UserParams::createForContiguous(rcppParams);
    std::function<hashing::ComplexHasher()> algorithmParams = [&userParams]() -> hashing::ComplexHasher {
        return createKMerComplexHasher(userParams.hashDim);
    };
    return countKMersSpecific<decltype(algorithmParams)>(
            sequences, alphabet, userParams, algorithmParams);
}

// [[Rcpp::export(".count_contiguous_kmers_string")]]
Rcpp::List count_contiguous_kmers_string(
        Rcpp::StringMatrix &sequenceMatrix,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countContiguousKMers(
            sequenceMatrix, alphabet, rcppParams
    );
}

// [[Rcpp::export(".count_contiguous_kmers_integer")]]
Rcpp::List count_contiguous_kmers_integer(
        Rcpp::IntegerMatrix &sequenceMatrix,
        Rcpp::IntegerVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countContiguousKMers(
            sequenceMatrix, alphabet, rcppParams
    );
}

// [[Rcpp::export(".count_contiguous_kmers_numeric")]]
Rcpp::List count_contiguous_kmers_numeric(
        Rcpp::NumericMatrix &sequenceMatrix,
        Rcpp::NumericVector &alphabet,
        Rcpp::Environment &rcppParams) {
    return countContiguousKMers(
            sequenceMatrix, alphabet, rcppParams
    );
}

// [[Rcpp::export(".count_contiguous_kmers_list")]]
Rcpp::List count_contiguous_kmers_list(
        Rcpp::List &sq,
        Rcpp::StringVector &alphabet,
        Rcpp::Environment& rcppParams) {
    return countContiguousKMers(
            sq, alphabet, rcppParams);
}
