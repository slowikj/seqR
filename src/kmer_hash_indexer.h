#ifndef KMER_HASH_INDEXER_H
#define KMER_HASH_INDEXER_H

#include "dictionary.h"
#include "kmer_counts_manager.h"
#include "hash/custom_hashers.h"
#include <vector>
#include <tuple>
#include <memory>

class KMerPositionInfo {
public:
    int seqNum;

    int position;

    KMerPositionInfo(int seqNum, int position) :
            seqNum(seqNum), position(position) {
    }

    KMerPositionInfo(const KMerPositionInfo &) = default;

    KMerPositionInfo &operator=(const KMerPositionInfo &) = default;

};

inline
std::tuple<Dictionary<std::vector<int>, int, vector_int_hasher>,
        std::vector<KMerPositionInfo>>
indexKMerHashes(const std::vector<KMerCountsManager> &kmerCounts) {
    Dictionary<std::vector<int>, int, vector_int_hasher> hashIndexer;
    std::vector<KMerPositionInfo> uniqueKMers;
    int currentIndex = 0;
    for (int seqNum = 0; seqNum < kmerCounts.size(); ++seqNum) {
        for (const auto &hashPair: kmerCounts[seqNum].getDictionary()) {
            auto hash = hashPair.first;
            auto seqStartPosition = hashPair.second.seqStartPosition;
            if (!hashIndexer.isPresent(hash)) {
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

#endif
