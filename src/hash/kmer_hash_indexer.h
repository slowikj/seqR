#ifndef KMER_HASH_INDEXER_H
#define KMER_HASH_INDEXER_H

#include "../dictionary.h"
#include <vector>
#include <tuple>
#include <memory>

struct KMerPositionInfo {
public:
  int seqNum;
  
  int position;
  
  KMerPositionInfo(int seqNum, int position):
    seqNum(seqNum), position(position) {
  }
  
  KMerPositionInfo(const KMerPositionInfo&) = default;
  
  KMerPositionInfo& operator=(const KMerPositionInfo&) = default;
  
};

std::tuple<Dictionary<std::vector<int>, int, vector_int_hasher>,
           std::vector<KMerPositionInfo>>
indexHashes(const std::vector<KMerCountsManager> kmerCounts);

#endif