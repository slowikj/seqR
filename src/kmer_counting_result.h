#pragma once

#include <Rcpp.h>

#include <string>
#include <tuple>
#include <vector>

#include "hash/custom_vector_hasher.h"
#include "hash/globals.h"

template <template <class K, class V, class...> class kmer_dictionary_t>
class KMerCountingResult {
 public:
  using hash_t = hashing::config::multidim_hash_t;
  using hasher_t = hashing::config::multidim_hasher_t;

  std::vector<std::string> kMerStrings;

  inline void increaseProcessedSequencesNum(int cnt) {
    this->processedSequencesNum += cnt;
  }

  inline bool addKMer(const hashing::config::multidim_hash_t &kMerHash,
                      int sequenceNum,
                      int count) {
    if (kMerHash2ColumnIndex.isPresent(kMerHash)) {
      addKMerCountsInfo(kMerHash2ColumnIndex[kMerHash], sequenceNum, count);
      return false;
    } else {
      int kMerIndex = kMerHash2ColumnIndex.size();
      kMerHash2ColumnIndex[kMerHash] = kMerIndex;
      addKMerCountsInfo(kMerIndex, sequenceNum, count);
      return true;
    }
  }

  inline Rcpp::List toRcppList() const {
    Rcpp::IntegerVector rcppSequenceNums = Rcpp::wrap(this->sequenceNums);
    Rcpp::IntegerVector rcppKMerIndices = Rcpp::wrap(this->kMerIndices);
    Rcpp::IntegerVector rcppKmerCounts = Rcpp::wrap(this->kMerCounts);
    Rcpp::StringVector rcppKMerStrings = Rcpp::wrap(this->kMerStrings);
    return Rcpp::List::create(
        Rcpp::Named(PROXY_ROWS_NAME) = rcppSequenceNums + 1,
        Rcpp::Named(PROXY_COLUMNS_NAME) = rcppKMerIndices + 1,
        Rcpp::Named(PROXY_VALUES_NAME) = rcppKmerCounts,
        Rcpp::Named(PROXY_COLUMN_NAMES_NAME) = rcppKMerStrings,
        Rcpp::Named(PROXY_PROCESSED_SEQUENCES_NUM_NAME) = this->processedSequencesNum);
  }

 private:
  const Rcpp::String PROXY_ROWS_NAME = "i";
  const Rcpp::String PROXY_COLUMNS_NAME = "j";
  const Rcpp::String PROXY_VALUES_NAME = "v";
  const Rcpp::String PROXY_COLUMN_NAMES_NAME = "names";
  const Rcpp::String PROXY_PROCESSED_SEQUENCES_NUM_NAME = "seqNum";

  std::vector<int> sequenceNums;

  std::vector<int> kMerIndices;

  std::vector<int> kMerCounts;

  kmer_dictionary_t<hash_t, int, hasher_t> kMerHash2ColumnIndex;

  int processedSequencesNum = 0;

  inline void addKMerCountsInfo(int kMerIndex, int sequenceNum, int kMerCount) {
    sequenceNums.push_back(sequenceNum + processedSequencesNum);
    kMerIndices.push_back(kMerIndex);
    kMerCounts.push_back(kMerCount);
  }
};
