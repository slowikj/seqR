#ifndef KMER_COUNTS_MANAGER_H
#define KMER_COUNTS_MANAGER_H

#include "dictionary.h"
#include "hash/custom_hashers.h"
#include <vector>
#include <memory>

struct KMerHashInfo {
  
  int cnt;
  
  int seqStartPosition;
  
  KMerHashInfo(int seqStartPosition):
    cnt(0), seqStartPosition(seqStartPosition) {
  }
  
  KMerHashInfo() = default;

  KMerHashInfo& operator=(const KMerHashInfo&) = default;
  
  ~KMerHashInfo() = default;
};

class KMerCountsManager {
public:
  
  KMerCountsManager() = default;
  
  KMerCountsManager(const KMerCountsManager&) = default;
  
  KMerCountsManager(KMerCountsManager&&) noexcept = default;
  
  KMerCountsManager& operator=(const KMerCountsManager&) = default;
  
  KMerCountsManager& operator=(KMerCountsManager&&) noexcept = default;
  
  void add(std::vector<int>&& hash, int position) {
    if(!this->dictionary.isPresent(hash)) {
      this->dictionary[hash] = KMerHashInfo(position);
    }
    this->dictionary[std::move(hash)].cnt++; 
  }
  
  const Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher>& getDictionary() const {
    return this->dictionary;
  }
  
private:
  Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher> dictionary;
  
};

#endif // KMER_COUNTS_MANAGER_H
