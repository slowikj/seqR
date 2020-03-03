#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include "hash/rolling_window.h"
#include "kmer_counts_manager.h"

template<class input_vector_t, class input_elem_t, class encoded_elem_t>
inline void updateKMerCounts(RollingWindow<input_vector_t, input_elem_t, encoded_elem_t>& rollingWindow,
                             KMerCountsManager& kmerCountsManager,
                             bool isPositionalKMer) {
  kmerCountsManager.add(
    isPositionalKMer ? rollingWindow.getWindowedPositionedHashes()
                     : rollingWindow.getWindowedHashes(),
    rollingWindow.currentBeginIndex()
  );
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t>
inline void countKMersForContiguousSeq(int k,
                                       int begin,
                                       int end,
                                       RollingWindow<input_vector_t, input_elem_t, encoded_elem_t>& rollingWindow,
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



#endif //KMER_COUNTER_H