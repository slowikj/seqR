#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include "hash/rolling_window.h"
#include "kmer_counts_manager.h"

void countKMersForContiguousSeq(int k,
                                int begin,
                                int end,
                                RollingWindow& rollingWindow,
                                KMerCountsManager& kmerCountsManager,
                                bool isPositionalKMer) {
  rollingWindow.resetIndex(begin);
  for(int i = 0; i < k; ++i) {
    rollingWindow.append();
  }
  for(int beginPosition = begin; beginPosition <= end; ++beginPosition) {
    kmerCountsManager.add(
      isPositionalKMer ? rollingWindow.getWindowedPositionedHashes()
                       : rollingWindow.getWindowedHashes(),
      beginPosition
    );
    rollingWindow.moveWindowRight();
  }
}

#endif //KMER_COUNTER_H