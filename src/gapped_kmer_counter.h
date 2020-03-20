#ifndef GAPPED_KMER_COUNTER_H
#define GAPPED_KMER_COUNTER_H

#include "Rcpp.h"
#include "alphabet_encoder.h"
#include "hash/polynomial_single_hasher.h"
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

std::vector<std::pair<int, int>> getContiguousIntervals(const Rcpp::IntegerVector& gaps);

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
  
  std::vector<int> getHashPositional(int begin, int end) const {
    std::vector<int> res = std::move(getHash(begin, end));
    for(int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      int P = polynomialHasherConfigs[hasherInd].P;
      int M = polynomialHasherConfigs[hasherInd].M;
      res[hasherInd] = static_cast<int>((static_cast<long long>(res[hasherInd]) * P + begin) % M);
    }
    return std::move(res);
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

#endif
