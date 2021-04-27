// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <vector>
#include "hash/globals.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_task_param_dispatcher.h"
#include "hash/complex_hasher.h"

inline hashing::ComplexHasher createKMerComplexHasher(int hashDim)
{
    std::vector<std::unique_ptr<hashing::SingleHasher>> singleHashers;
    for (int i = 0; i < hashDim; ++i)
    {
        singleHashers.emplace_back(new hashing::PolynomialSingleHasher(
            hashing::PolynomialSingleHasherConfig(
                hashing::config::hashPrimes[i].first,
                hashing::config::hashPrimes[i].second)));
    }
    hashing::ComplexHasher complexHasher(std::move(singleHashers));
    return complexHasher;
}

template <class sequences_t,
          class alphabet_t>
Rcpp::List countContiguousKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    Rcpp::Environment &rcppParams)
{
    auto userParams = UserParams::createForContiguous(rcppParams);
    std::function<hashing::ComplexHasher()> algorithmParams = [&userParams]() -> hashing::ComplexHasher {
        return createKMerComplexHasher(userParams.hashDim);
    };
    return countKMers<sequences_t, alphabet_t, decltype(algorithmParams)>(
        sequences, alphabet, userParams, algorithmParams);
}

// [[Rcpp::export(".count_contiguous_kmers_string_vector")]]
Rcpp::List count_contiguous_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams)
{
    return countContiguousKMers(sq, alphabet, rcppParams);
}
