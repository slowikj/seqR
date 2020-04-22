// [[Rcpp::plugins("cpp17")]]
#include "kmer_hash_indexer.h"

std::tuple<Dictionary<std::vector<int>, int, vector_int_hasher>,
           std::vector<KMerPositionInfo>>
indexKMerHashes(const std::vector<KMerCountsManager>& kmerCounts) {
  Dictionary<std::vector<int>, int, vector_int_hasher> hashIndexer;
  std::vector<KMerPositionInfo> uniqueKMers;
  int currentIndex = 0;
  for(int seqNum = 0; seqNum < kmerCounts.size(); ++seqNum) {
    for(const auto& hashPair: kmerCounts[seqNum].getDictionary()) {
      auto hash = hashPair.first;
      auto seqStartPosition = hashPair.second.seqStartPosition;
      if(!hashIndexer.isPresent(hash)) {
        hashIndexer[hash] = currentIndex++;
        uniqueKMers.emplace_back(seqNum, seqStartPosition);
      }
    }
  }
  return {
    std::move(hashIndexer),
    std::move(uniqueKMers)
  };
}
