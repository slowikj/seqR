// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>

#include "hash/complex_hasher.h"
#include "hash/globals.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_task_param_dispatcher.h"
#include "merge_kmer_results.h"

// ------------- EXPORTED -------------

Rcpp::List count_contiguous_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams);

Rcpp::List count_contiguous_kmers_string_list(
    Rcpp::List &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams);

Rcpp::List count_gapped_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams);

Rcpp::List count_gapped_kmers_string_list(
    Rcpp::List &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams);

Rcpp::List merge_kmer_results(
    Rcpp::List resList);

// ------------- HELPERS: CONTIGUOUS K-MERS COUNTING -------------

inline hashing::ComplexHasher createKMerComplexHasher(
    int hashDim);

template <class sequences_t,
          class alphabet_t>
inline Rcpp::List countContiguousKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    Rcpp::Environment &rcppParams);

// ------------- HELPERS: GAPPED K-MERS COUNTING -------------

inline std::vector<hashing::PolynomialSingleHasherConfig> getGappedKMerHasherConfigs(
    int hashDim);

template <class sequences_t,
          class alphabet_t>
Rcpp::List countGappedKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    Rcpp::Environment &rcppParams);

// ------------- IMPLEMENTATION -------------

// ------------- CONTIGUOUS K-MER COUNTING -------------

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
          class alphabet_t>
Rcpp::List countContiguousKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    Rcpp::Environment &rcppParams) {
  auto userParams = UserParams::createForContiguous(rcppParams);
  std::function<hashing::ComplexHasher()> algorithmParams = [&userParams]() -> hashing::ComplexHasher {
    return createKMerComplexHasher(userParams.hashDim);
  };
  return countKMers<sequences_t, alphabet_t, decltype(algorithmParams)>(
      sequences, alphabet, userParams, algorithmParams);
}

// [[Rcpp::export(".cpp_count_contiguous_kmers_string_vector")]]
Rcpp::List count_contiguous_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams) {
  return countContiguousKMers(sq, alphabet, rcppParams);
}

// [[Rcpp::export(".cpp_count_contiguous_kmers_string_list")]]
Rcpp::List count_contiguous_kmers_string_list(
    Rcpp::List &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams) {
  return countContiguousKMers(sq, alphabet, rcppParams);
}

// ------------- GAPPED K-MER COUNTING -------------

inline std::vector<hashing::PolynomialSingleHasherConfig> getGappedKMerHasherConfigs(int hashDim) {
  std::vector<hashing::PolynomialSingleHasherConfig> res;
  for (int i = 0; i < hashDim; ++i) {
    res.emplace_back(hashing::config::hashPrimes[i].first, hashing::config::hashPrimes[i].second);
  }
  return res;
}

template <class sequences_t,
          class alphabet_t>
Rcpp::List countGappedKMers(
    sequences_t &sequences,
    alphabet_t &alphabet,
    Rcpp::Environment &rcppParams) {
  auto userParams = UserParams::createForGapped(rcppParams);
  auto hasherConfigs = getGappedKMerHasherConfigs(userParams.hashDim);
  return countKMers<sequences_t, alphabet_t, decltype(hasherConfigs)>(
      sequences, alphabet, userParams, hasherConfigs);
}

// [[Rcpp::export(".cpp_count_gapped_kmers_string_vector")]]
Rcpp::List count_gapped_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams) {
  return countGappedKMers(sq, alphabet, rcppParams);
}

// [[Rcpp::export(".cpp_count_gapped_kmers_string_list")]]
Rcpp::List count_gapped_kmers_string_list(
    Rcpp::List &sq,
    Rcpp::StringVector &alphabet,
    Rcpp::Environment &rcppParams) {
  return countGappedKMers(sq, alphabet, rcppParams);
}

// ------------- MERGING K-MER RESULTS -------------

// [[Rcpp::export(".cpp_merge_kmer_results")]]
Rcpp::List merge_kmer_results(Rcpp::List resList) {
  return resultsMerging::mergeKMerResults(resList);
}
