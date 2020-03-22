#ifndef GAPPED_KMER_COUNTER_H
#define GAPPED_KMER_COUNTER_H

#include "Rcpp.h"
#include "alphabet_encoder.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_counting_common.h"
#include "kmer_counts_manager.h"
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

std::vector<std::pair<int, int>> getContiguousIntervals(const Rcpp::IntegerVector& gaps);

std::size_t getTotalKMerSize(const Rcpp::IntegerVector& gaps);

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<int> prepareNotAllowedItemsPrefixCount(
    input_vector_t& sequence,
    const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding) {
  std::vector<int> res;
  res.reserve(sequence.size());
  for(int i = 0; i < sequence.size(); ++i) {
    bool isNotPresent = !alphabetEncoding.isAllowed(sequence[i]);
    res[i] = (i == 0) ?
      isNotPresent :
      res[i - 1] + isNotPresent;
  }
  return std::move(res);
}

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
class PrefixSequencePolynomialHasher {
public:
  PrefixSequencePolynomialHasher(input_vector_t& sequence,
                                 const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
                                 std::vector<PolynomialSingleHasherConfig>&& polynomialHasherConfigs):
    polynomialHasherConfigs(std::move(polynomialHasherConfigs)) {
    computePrefixValues(sequence, alphabetEncoding);
  }
  
  std::vector<int> getHash(int begin, int end) const {
    std::vector<int> res(polynomialHasherConfigs.size());
    for(int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      int M = polynomialHasherConfigs[hasherInd].M;
      res[hasherInd] = static_cast<int>((prefixComplexHashes[end + 1][hasherInd]
        - static_cast<long long>(prefixComplexHashes[begin][hasherInd]) * prefixP[end - begin][hasherInd] + M) % M);
    }
    return std::move(res);
  }
  
  int getHasherP(int hasherIndex, int power = 1) const {
    return this->prefixP[power][hasherIndex];
  }
  
  int getHasherM(int hasherIndex) const {
    return this->polynomialHasherConfigs[hasherIndex].M;
  }

  int getHashersNum() const {
    return this->polynomialHasherConfigs.size();
  }
  
private:
  std::vector<PolynomialSingleHasherConfig> polynomialHasherConfigs;
  std::vector<std::vector<int>> prefixP;
  std::vector<std::vector<int>> prefixComplexHashes;
  
  void computePrefixValues(input_vector_t& sequence,
                           const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding) {
    initPrefixP(sequence.size(), polynomialHasherConfigs.size());
    initPrefixComplexHashes(sequence.size(), polynomialHasherConfigs.size());
    for(int seqInd  = 0; seqInd < sequence.size(); ++seqInd) {
      appendPrefixValues(sequence, alphabetEncoding, seqInd);
    }
  }
  
  void initPrefixP(int sequenceLength, int hashersNum) {
    prefixP.reserve(sequenceLength);
    prefixP.push_back(std::move(std::vector<int>(hashersNum, 1)));
  }
  
  void initPrefixComplexHashes(int sequenceLength, int hashersNum) {
    prefixComplexHashes.reserve(sequenceLength);
    prefixComplexHashes.push_back(
      std::move(std::vector<int>(hashersNum))
    );
  }
  
  void appendPrefixValues(input_vector_t& sequence,
                          const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
                          int seqInd) {
    encoded_elem_t encodedElem = alphabetEncoding.isAllowed(sequence[seqInd]) ?
      alphabetEncoding.encode(sequence[seqInd]) :
      alphabetEncoding.getNotAllowedEncodingNum();
    
    appendCurrentComplexHash();
    appendCurrentPowerP();
  }
  
  void appendCurrentComplexHash() {
    std::vector<int> prefixHash;
    std::transform(std::begin(polynomialHasherConfigs), std::end(polynomialHasherConfigs),
                   std::back_inserter(prefixHash),
                   [this](const PolynomialSingleHasherConfig& hasherConfig) -> int {
                     return static_cast<int>(
                       (static_cast<long long>(prefixComplexHashes.back()) * hasherConfig.P) %  hasherConfig.M);
                   }
    );
    prefixComplexHashes.push_back(std::move(prefixHash));
  }
  
  void appendCurrentPowerP() {
    std::vector<int> powersP;
    std::transform(std::begin(polynomialHasherConfigs), std::end(polynomialHasherConfigs),
                   std::back_inserter(powersP),
                   [this](const PolynomialSingleHasherConfig& hasherConfig) -> int {
                     return static_cast<int>(
                       (static_cast<long long>(prefixP.back()) * hasherConfig.P) % hasherConfig.M);
                   }
    );
    prefixP.push_back(std::move(powersP));
  }
};

bool isGappedKMerAllowed(const std::vector<std::pair<int,int>>& contiguousKMerIntervals,
                         const std::vector<int>& notAllowedItemsPrefixCount);

int getIntervalLength(const std::pair<int, int>& interval);

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<int> getGappedKMerHashNotPositional(
    int beginPosition,
    const PrefixSequencePolynomialHasher<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>& seqHasher,
    const std::vector<std::pair<int,int>>& contiguousKMerIntervals) {
  std::vector<int> res(seqHasher.getHashersNum());
  for(const auto& kmerInterval: contiguousKMerIntervals) {
    auto intervalHash = std::move(seqHasher.getHash(
      kmerInterval.first + beginPosition,
      kmerInterval.second + beginPosition));
    int intervalLength = getIntervalLength(kmerInterval);
    for(int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      int powerP = seqHasher.getHasherP(hasherInd, intervalLength);
      int M = seqHasher.getHasherM(hasherInd);
      res[hasherInd] = static_cast<int>((
        static_cast<long long>(res[hasherInd]) * powerP + intervalHash[hasherInd]) % M);
    }
  }
  return std::move(res);
}

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<int> getGappedKMerHash(
  int beginPosition,
  const PrefixSequencePolynomialHasher<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>& seqHasher,
  const std::vector<std::pair<int,int>>& contiguousKMerIntervals,
  bool isPositionalKMer) {
  std::vector<int> res = std::move(
    getGappedKMerHashNotPositional<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
     beginPosition, seqHasher, contiguousKMerIntervals 
  ));
  if(isPositionalKMer) {
    for(int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      res[hasherInd] = static_cast<int>(
        (static_cast<long long>(res[hasherInd]) * seqHasher.getHasherP(hasherInd)
           + beginPosition) % seqHasher.getHasherM(hasherInd));
    }
  }
  return std::move(res);
}

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
KMerCountsManager countGappedKMers(const Rcpp::IntegerVector& gaps,
                                   std::size_t totalKMerSize,
                                   input_vector_t& sequence,
                                   AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
                                   bool isPositionalKMer,
                                   std::vector<PolynomialSingleHasherConfig>&& hasherConfigs) {
  std::vector<std::pair<int,int>> contiguousIntervals = getContiguousIntervals(gaps);
  PrefixSequencePolynomialHasher<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t> sequenceHasher(
      sequence, alphabetEncoding, std::move(hasherConfigs)
  );
  
  std::vector<int> notAllowedItemsPrefixCount = std::move(
    prepareNotAllowedItemsPrefixCount<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
      sequence,
      alphabetEncoding
    )
  );

  KMerCountsManager kmerCountsManager;
  for(int seqInd = 0; seqInd < sequence.size() - totalKMerSize + 1; ++seqInd) {
    if(isGappedKMerAllowed(contiguousIntervals, notAllowedItemsPrefixCount)) {
      auto hash = std::move(getGappedKMerHash(seqInd, sequenceHasher, contiguousIntervals, isPositionalKMer));
      kmerCountsManager.add(hash, seqInd);
    }
  }
  
  return std::move(kmerCountsManager);
}

template<class input_matrix_t, class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<KMerCountsManager> parallelComputeGappedKMersCounts(
  const Rcpp::IntegerVector& gaps,
  bool isPositionalKMer,
  input_matrix_t& sequenceMatrix,
  AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding) {
  std::size_t totalKMerSize = getTotalKMerSize(gaps);
  return std::move(parallelComputeKMerCounts<input_matrix_t, input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
    sequenceMatrix,
    [&gapsVector = std::as_const(gaps), isPositionalKMer, &alphabetEncoding, &totalKMerSize]
    (input_vector_t& v) -> KMerCountsManager {
      return countGappedKMers<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
          gapsVector,
          totalKMerSize,
          v,
          alphabetEncoding,
          isPositionalKMer
      );
    }
  ));
}

#endif
