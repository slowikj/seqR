#pragma once

#include <memory>
#include <string>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>

#include <algorithm>
#include <functional>

#include "kmer_task_config.h"
#include "utils.h"

namespace stringsCreator {

class KMerPositionInfo;

template <class encoded_sequence_t>
class KMerStringCreatorForSequence;

template <class encoded_sequences_list_t>
class KMerStringsCreatorWorker;

template <class encoded_sequences_list_t>
inline void generate(
    const std::vector<KMerPositionInfo> &indexedKMers,
    const KMerTaskConfig<encoded_sequences_list_t> &kMerTaskConfig,
    std::vector<std::string> &resultStrings);

// ------------------ IMPLEMENTATION ------------------

template <class encoded_sequences_list_t>
inline void generate(
    const std::vector<KMerPositionInfo> &indexedKMers,
    const KMerTaskConfig<encoded_sequences_list_t> &kMerTaskConfig,
    std::vector<std::string> &resultStrings) {
  KMerStringsCreatorWorker<encoded_sequences_list_t> worker(indexedKMers, kMerTaskConfig, resultStrings);
  RcppParallel::parallelFor(0, indexedKMers.size(), worker);
}

class KMerPositionInfo {
 public:
  int seqNum;

  int position;

  KMerPositionInfo(int seqNum, int position) : seqNum(seqNum), position(position) {
  }

  KMerPositionInfo(const KMerPositionInfo &) = default;

  KMerPositionInfo &operator=(const KMerPositionInfo &) = default;
};

template <class encoded_sequences_list_t>
class KMerStringsCreatorWorker : public RcppParallel::Worker {
 public:
  KMerStringsCreatorWorker(
      const std::vector<KMerPositionInfo> &kMersToGenerate,
      const KMerTaskConfig<encoded_sequences_list_t> &kMerTaskConfig,
      std::vector<std::string> &resultStrings)
      : kMersToGenerate(kMersToGenerate),
        kMerTaskConfig(kMerTaskConfig),
        gapsAccumulated(util::getGapsAccumulated(kMerTaskConfig.userParams.gaps)),
        resultStrings(resultStrings),
        resultOffset(resultStrings.size()) {
    resultStrings.resize(resultOffset + kMersToGenerate.size());
    prepareKMerStringsCreators();
    prepareCreateKMerFunc();
  }

  inline void operator()(std::size_t begin, std::size_t end) override {
    for (int i = begin; i < end; ++i) {
      this->resultStrings[i + resultOffset] = createKMerFunc(
          this->kMersToGenerate[i].seqNum, this->kMersToGenerate[i].position);
    }
  }

 private:
  const std::vector<KMerPositionInfo> &kMersToGenerate;
  const KMerTaskConfig<encoded_sequences_list_t> &kMerTaskConfig;
  std::vector<KMerStringCreatorForSequence<typename encoded_sequences_list_t::Entry>> kmerStringCreators;
  std::function<std::string(int, int)> createKMerFunc;
  std::vector<int> gapsAccumulated;
  std::vector<std::string> &resultStrings;
  int resultOffset;

  inline void prepareKMerStringsCreators() {
    std::size_t sequencesNum = kMerTaskConfig.encodedSequencesList.size();
    this->kmerStringCreators.reserve(sequencesNum);
    for (int i = 0; i < sequencesNum; ++i) {
      auto seq = kMerTaskConfig.encodedSequencesList[i];
      this->kmerStringCreators.emplace_back(
          std::move(seq),
          kMerTaskConfig.userParams.gaps, gapsAccumulated,
          kMerTaskConfig.kMerItemSeparator, kMerTaskConfig.kMerSectionSeparator);
    }
  }

  inline void prepareCreateKMerFunc() {
    if (kMerTaskConfig.userParams.positional) {
      this->createKMerFunc = [this](int seqNum, int pos) -> std::string {
        return kmerStringCreators[seqNum].getPositional(pos);
      };
    } else {
      this->createKMerFunc = [this](int seqNum, int pos) -> std::string {
        return kmerStringCreators[seqNum].get(pos);
      };
    }
  }
};

template <class encoded_sequence_t>
class KMerStringCreatorForSequence {
 public:
  KMerStringCreatorForSequence(
      encoded_sequence_t &&sequence,
      const std::vector<int> &gaps,
      const std::vector<int> &gapsAccumulated,
      const std::string &itemSeparator,
      const std::string &sectionSeparator)
      : sequence(std::move(sequence)),
        itemSeparator(itemSeparator),
        sectionSeparator(sectionSeparator),
        kmerInfoSuffix(prepareKMerInfoSuffix(gaps)),
        gapsAccumulated(gapsAccumulated) {
  }

  inline std::string get(std::size_t begin) const {
    int totalSize = getTotalSize(begin, itemSeparator.size());
    std::string res;
    res.reserve(totalSize);
    res += sequence.decode(begin);
    for (const int &accGap : this->gapsAccumulated) {
      res += itemSeparator + sequence.decode(begin + accGap);
    }
    return kmerInfoSuffix.size() > 0
               ? res + sectionSeparator + kmerInfoSuffix
               : res;
  }

  inline std::string getPositional(std::size_t begin) const {
    std::string withoutPositionString = this->get(begin);
    return std::to_string(begin + 1) + sectionSeparator + withoutPositionString;
  }

 private:
  encoded_sequence_t sequence;
  std::string itemSeparator;
  std::string sectionSeparator;
  std::string kmerInfoSuffix;
  const std::vector<int> &gapsAccumulated;

  inline std::size_t getTotalSize(std::size_t begin, int separatorLength) const {
    return std::accumulate(
        std::begin(this->gapsAccumulated),
        std::end(this->gapsAccumulated),
        sequence.decode(begin).size(),
        [this, &begin, &separatorLength](int r, int accGap) {
          return r + sequence.decode(begin + accGap).size() + separatorLength;
        });
  }

  inline std::string prepareKMerInfoSuffix(const std::vector<int> &gaps) const {
    if (gaps.empty()) {
      return "";
    }
    std::string res;
    int approximateResSize = gaps.size() + (gaps.size() - 1) * itemSeparator.size();
    res.reserve(approximateResSize);
    res += std::to_string(gaps[0]);
    for (int gaps_i = 1; gaps_i < gaps.size(); ++gaps_i) {
      res += itemSeparator;
      res += std::to_string(gaps[gaps_i]);
    }
    return res;
  }
};

}  // namespace stringsCreator
