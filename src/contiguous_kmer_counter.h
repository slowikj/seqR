#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "hash/complex_hasher.h"
#include "hash/rolling_window.h"

namespace contiguousKMer {

template <class encoded_sequence_t,
          class kmer_manager_t>
inline kmer_manager_t count(
    const encoded_sequence_t &sequence,
    int k,
    bool isPositionalKMer,
    bool withKMerCounts,
    hashing::ComplexHasher &&complexHasher);

template <class encoded_sequence_t>
inline std::vector<int> computeNotAllowedPositions(const encoded_sequence_t &sequence);

template <class encoded_sequence_t,
          class kmer_manager_t>
inline void countKMersForContiguousSeq(
    hashing::RollingWindow<encoded_sequence_t> &rollingWindow,
    kmer_manager_t &kMerManager,
    int k,
    int begin,
    int end,
    bool isPositionalKMer);

template <class encoded_sequence_t,
          class kmer_manager_t>
inline void updateKMers(
    hashing::RollingWindow<encoded_sequence_t> &rollingWindow,
    kmer_manager_t &kMerManager,
    bool isPositionalKMer);

// ------------------ IMPLEMENTATION ------------------

template <class encoded_sequence_t,
          class kmer_manager_t>
inline kmer_manager_t count(
    const encoded_sequence_t &sequence,
    std::size_t k,
    bool isPositionalKMer,
    hashing::ComplexHasher &&complexHasher) {
  kmer_manager_t kMerManager;
  hashing::RollingWindow<encoded_sequence_t> rollingWindow(
      sequence,
      std::move(complexHasher));
  auto notAllowedSequencePositions = computeNotAllowedPositions(sequence);
  for (std::size_t i = 0; i < notAllowedSequencePositions.size() - 1; ++i) {
    std::size_t allowedItemsBetween = notAllowedSequencePositions[i + 1] - notAllowedSequencePositions[i] - 1;
    if (allowedItemsBetween >= k) {
      std::size_t begin = notAllowedSequencePositions[i] + 1;
      std::size_t end = notAllowedSequencePositions[i + 1] - 1;
      countKMersForContiguousSeq<encoded_sequence_t, kmer_manager_t>(
          rollingWindow, kMerManager, k, begin, end, isPositionalKMer);
    }
  }
  return kMerManager;
}

template <class encoded_sequence_t>
inline std::vector<int> computeNotAllowedPositions(
    const encoded_sequence_t &sequence) {
  int leftSentinel = -1;
  int rightSentinel = sequence.size();

  if (sequence.areAllElementsAllowed()) {
    return {leftSentinel, rightSentinel};
  }

  std::vector<int> res;
  res.push_back(leftSentinel);
  for (std::size_t seq_i = 0; seq_i < sequence.size(); ++seq_i) {
    if (!sequence.isAllowed(seq_i)) {
      res.push_back(seq_i);
    }
  }
  res.push_back(rightSentinel);
  return res;
}

template <class encoded_sequence_t,
          class kmer_manager_t>
inline void countKMersForContiguousSeq(
    hashing::RollingWindow<encoded_sequence_t> &rollingWindow,
    kmer_manager_t &kMerManager,
    int k,
    int begin,
    int end,
    bool isPositionalKMer) {
  rollingWindow.resetIndex(begin);
  for (int i = 0; i < k; ++i) {
    rollingWindow.append();
  }
  for (int beginPosition = begin; beginPosition < end - k + 1; ++beginPosition) {
    updateKMers(rollingWindow, kMerManager, isPositionalKMer);
    rollingWindow.moveWindowRight();
  }
  updateKMers(rollingWindow, kMerManager, isPositionalKMer);
}

template <class encoded_sequence_t,
          class kmer_manager_t>
inline void updateKMers(
    hashing::RollingWindow<encoded_sequence_t> &rollingWindow,
    kmer_manager_t &kMerManager,
    bool isPositionalKMer) {
  kMerManager.add(
      isPositionalKMer ? rollingWindow.getWindowedPositionedHashes()
                       : rollingWindow.getWindowedHashes(),
      rollingWindow.currentBeginIndex());
}

}  // namespace contiguousKMer
