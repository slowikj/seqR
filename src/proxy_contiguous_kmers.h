#pragma once

// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>

#include "hash/complex_hasher.h"
#include "hash/globals.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_task_param_dispatcher.h"

inline hashing::ComplexHasher createKMerComplexHasher(int hashDim) {
  std::vector<std::unique_ptr<hashing::SingleHasher>> singleHashers;
  for (int i = 0; i < hashDim; ++i) {
    singleHashers.emplace_back(new hashing::PolynomialSingleHasher(
        hashing::PolynomialSingleHasherConfig(
            hashing::config::hashPrimes[i].first,
            hashing::config::hashPrimes[i].second)));
  }
  hashing::ComplexHasher complexHasher(std::move(singleHashers));
  return complexHasher;
}

template <class sequences_t,
          class kmer_alphabet_t>
inline Rcpp::List countContiguousKMers(
    sequences_t &sequences,
    kmer_alphabet_t &kmerAlphabet,
    Rcpp::Environment &rcppParams) {
  auto userParams = UserParams::createForContiguous(rcppParams);
  std::function<hashing::ComplexHasher()> algorithmParams = [&userParams]() -> hashing::ComplexHasher {
    return createKMerComplexHasher(userParams.hashDim);
  };
  return countKMers<sequences_t, kmer_alphabet_t, decltype(algorithmParams)>(
      sequences, kmerAlphabet, userParams, algorithmParams);
}
