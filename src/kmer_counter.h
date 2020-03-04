#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include "hash/rolling_window.h"
#include "hash/complex_hasher.h"
#include "hash/single_hasher.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_counts_manager.h"
#include <vector>

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
inline void updateKMerCounts(RollingWindow<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>& rollingWindow,
                             KMerCountsManager& kmerCountsManager,
                             bool isPositionalKMer) {
  kmerCountsManager.add(
    isPositionalKMer ? rollingWindow.getWindowedPositionedHashes()
                     : rollingWindow.getWindowedHashes(),
    rollingWindow.currentBeginIndex()
  );
}

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
inline void countKMersForContiguousSeq(int k,
                                       int begin,
                                       int end,
                                       RollingWindow<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>& rollingWindow,
                                       KMerCountsManager& kmerCountsManager,
                                       bool isPositionalKMer) {
  rollingWindow.resetIndex(begin);
  for(int i = 0; i < k; ++i) {
    rollingWindow.append();
  }
  for(int beginPosition = begin; beginPosition < end; ++beginPosition) {
    updateKMerCounts(rollingWindow, kmerCountsManager, isPositionalKMer);
    rollingWindow.moveWindowRight();
  }
  updateKMerCounts(rollingWindow, kmerCountsManager, isPositionalKMer);
}

inline ComplexHasher createComplexHasher() {
  std::vector<std::unique_ptr<SingleHasher>> singleHashers;
  singleHashers.emplace_back(new PolynomialSingleHasher(101, 1e9 + 33));
  singleHashers.emplace_back(new PolynomialSingleHasher(97, 1e9 + 7));
  ComplexHasher complexHasher(std::move(singleHashers));
  return complexHasher;
}

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
inline std::vector<int> computeNotAllowedPositions(
    AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
    input_vector_t& sequence) {
  std::vector<int> res;
  res.push_back(-1); // left sentinel
  for(int seq_i = 0; seq_i < sequence.size(); ++seq_i) {
    if(!alphabetEncoding.isAllowed(sequence[seq_i])) {
      res.push_back(seq_i);
    }
  }
  res.push_back(sequence.size()); // right sentinel
  return res;
}

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
inline KMerCountsManager countKMers(int k,
                                    input_vector_t& sequence,
                                    AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
                                    bool isPositionalKMer) {
  KMerCountsManager kmerCountsManager;
  RollingWindow<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t> rollingWindow(
      sequence, std::move(createComplexHasher()), alphabetEncoding
  );
  auto notAllowedSequencePositions = computeNotAllowedPositions(alphabetEncoding, sequence);
  for(int i = 0; i < notAllowedSequencePositions.size() - 1; ++i) {
    int allowedItemsBetween = notAllowedSequencePositions[i + 1] - notAllowedSequencePositions[i];
    if(allowedItemsBetween >= k) {
      int begin = notAllowedSequencePositions[i] + 1;
      int end = notAllowedSequencePositions[i + 1] - 1;
      countKMersForContiguousSeq<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
          k, begin, end, rollingWindow, kmerCountsManager, isPositionalKMer
      );
    }
  }
  return kmerCountsManager;
}

#endif //KMER_COUNTER_H
