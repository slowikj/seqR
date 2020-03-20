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
                                 std::vector<PolynomialSingleHasher>&& polynomialHashers):
    polynomialHashers(std::move(polynomialHashers)) {
    computePrefixValues(sequence, alphabetEncoding);
  }
  
  std::vector<int> getHash(int begin, int end) const {
    std::vector<int> res(polynomialHashers.size());
    for(int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      int M = polynomialHashers[hasherInd].getM();
      res[hasherInd] = static_cast<int>((prefixComplexHashes[end + 1][hasherInd]
        - static_cast<long long>(prefixComplexHashes[begin][hasherInd]) * prefixP[end - begin][hasherInd] + M) % M);
    }
    return std::move(res);
  }
  
  std::vector<int> getHashPositional(int begin, int end) const {
    std::vector<int> res = std::move(getHash(begin, end));
    for(int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      int P = polynomialHashers[hasherInd].getP();
      int M = polynomialHashers[hasherInd].getM();
      res[hasherInd] = static_cast<int>((static_cast<long long>(res[hasherInd]) * P + begin) % M);
    }
    return std::move(res);
  }
  
private:
  std::vector<PolynomialSingleHasher> polynomialHashers;
  std::vector<std::vector<int>> prefixP;
  std::vector<std::vector<int>> prefixComplexHashes;
  
  void computePrefixValues(input_vector_t& sequence,
                           const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding) {
    prefixP.reserve(sequence.size());
    prefixComplexHashes.reserve(sequence.size() + 1);
    prefixComplexHashes.push_back(std::move(std::vector<int>(polynomialHashers.size())));
    for(int seqInd  = 0; seqInd < sequence.size(); ++seqInd) {
      appendPrefixValues(sequence, alphabetEncoding, seqInd);
    }
  }
  
  void appendPrefixValues(input_vector_t& sequence,
                          const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
                          int seqInd) {
    encoded_elem_t encodedElem = alphabetEncoding.isAllowed(sequence[seqInd]) ?
      alphabetEncoding.encode(sequence[seqInd]) :
      alphabetEncoding.getNotAllowedEncodingNum();
    
    hashElement(encodedElem);
    appendCurrentComplexHash();
    appendCurrentPowerP();
  }
  
  void hashElement(const encoded_elem_t& encodedElem) {
    std::for_each(std::begin(polynomialHashers), std::end(polynomialHashers),
                  [&encodedElem](PolynomialSingleHasher& hasher) {
                    hasher.append(encodedElem);
                  }
    );
  }
  
  void appendCurrentComplexHash() {
    std::vector<int> prefixHash;
    std::transform(std::begin(polynomialHashers), std::end(polynomialHashers),
                   std::back_inserter(prefixHash),
                   [](const PolynomialSingleHasher& hasher) -> std::vector<int> {
                     return hasher.getHash();
                   }
    );
    prefixComplexHashes.push_back(std::move(prefixHash));
  }
  
  void appendCurrentPowerP() {
    std::vector<int> powersP;
    std::transform(std::begin(polynomialHashers), std::end(polynomialHashers),
                   std::back_inserter(powersP),
                   [](const PolynomialSingleHasher& hasher) -> std::vector<int> {
                     return hasher.getCurrentPowerP();
                   }
    );
    prefixP.push_back(std::move(powersP));
  }
};

#endif
