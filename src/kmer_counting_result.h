#ifndef SEQR_KMER_COUNTING_RESULT_H
#define SEQR_KMER_COUNTING_RESULT_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>
#include "hash/custom_hashers.h"

class KMerCountingResult {
public:

    const Rcpp::String PROXY_ROWS_NAME = "i";
    const Rcpp::String PROXY_COLUMNS_NAME = "j";
    const Rcpp::String PROXY_VALUES_NAME = "v";
    const Rcpp::String PROXY_COLUMN_NAMES_NAME = "names";
    const Rcpp::String PROXY_PROCESSED_SEQUENCES_NAME = "processed_sequences";

    std::vector<int> sequenceNums;

    std::vector<int> kMerIndices;

    std::vector<int> kMerCounts;

    std::vector<std::string> kMerStrings;

    int processedSequences = 0;

    inline bool addKMer(const std::vector<int> &kMerHash,
                        int sequenceNum,
                        int count) {
        if (kMerHash2ColumnIndex.find(kMerHash) != kMerHash2ColumnIndex.end()) {
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
                Rcpp::Named(PROXY_PROCESSED_SEQUENCES_NAME) = this->processedSequences);

    }

private:
    std::unordered_map<std::vector<int>, int> kMerHash2ColumnIndex;

    inline void addKMerCountsInfo(int kMerIndex, int sequenceNum, int kMerCount) {
        sequenceNums.push_back(sequenceNum + processedSequences);
        kMerIndices.push_back(kMerIndex);
        kMerCounts.push_back(kMerCount);
    }

};

#endif //SEQR_KMER_COUNTING_RESULT_H