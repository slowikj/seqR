#pragma once

// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>

#include "hash/complex_hasher.h"
#include "hash/globals.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_task_param_dispatcher.h"

inline std::vector<hashing::PolynomialSingleHasherConfig> getGappedKMerHasherConfigs(int hashDim) {
  std::vector<hashing::PolynomialSingleHasherConfig> res;
  for (int i = 0; i < hashDim; ++i) {
    res.emplace_back(hashing::config::hashPrimes[i].first, hashing::config::hashPrimes[i].second);
  }
  return res;
}

template <class sequences_t,
          class kmer_alphabet_t>
inline Rcpp::List countGappedKMers(
    sequences_t &sequences,
    kmer_alphabet_t &kmerAlphabet,
    Rcpp::Environment &rcppParams) {
  auto userParams = UserParams::createForGapped(rcppParams);
  auto hasherConfigs = getGappedKMerHasherConfigs(userParams.hashDim);
  return countKMers<sequences_t, kmer_alphabet_t, decltype(hasherConfigs)>(
      sequences, kmerAlphabet, userParams, hasherConfigs);
}
