#ifndef KMER_STRINGS_CREATOR_H
#define KMER_STRINGS_CREATOR_H

#include<string>
#include<memory>
#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include <functional>
#include <algorithm>
#include "hash/complex_hasher.h"
#include "utils.h"
#include "kmer_task_config.h"

namespace stringsCreator {

    class KMerPositionInfo;

    template<class input_vector_t, class input_elem_t>
    class KMerStringCreatorForSequence;

    template<class input_vector_t, class input_elem_t>
    class KMerStringsCreatorWorker;

    template<class input_vector_t, class input_elem_t>
    inline void generate(
            const std::vector<KMerPositionInfo> &indexedKMers,
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            std::vector<std::string> &resultStrings);

    template<class input_vector_t, class input_elem_t>
    inline void generate(
            const std::vector<KMerPositionInfo> &indexedKMers,
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            std::vector<std::string> &resultStrings) {
        KMerStringsCreatorWorker<input_vector_t, input_elem_t> worker(indexedKMers, kMerTaskConfig, resultStrings);
        if (kMerTaskConfig.parallelMode) {
            RcppParallel::parallelFor(0, indexedKMers.size(), worker);
        } else {
            worker(0, indexedKMers.size());
        }
    }

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

    template<class input_vector_t, class input_elem_t>
    class KMerStringsCreatorWorker : public RcppParallel::Worker {
    public:
        KMerStringsCreatorWorker(const std::vector<KMerPositionInfo> &kMersToGenerate,
                                 KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                                 std::vector<std::string> &resultStrings) :
                kMersToGenerate(kMersToGenerate),
                kMerTaskConfig(kMerTaskConfig),
                gapsAccumulated(std::move(util::getGapsAccumulated(kMerTaskConfig.gaps))),
                resultOffset(resultStrings.size()),
                resultStrings(resultStrings) {
            resultStrings.resize(resultOffset + kMersToGenerate.size());
            prepareKMerStringsCreators();
            prepareCreateKMerFunc();
        }

        inline void operator()(std::size_t begin, std::size_t end) override {
            for (int i = begin; i < end; ++i) {
                this->resultStrings[i + resultOffset] = std::move(
                        createKMerFunc(this->kMersToGenerate[i].seqNum, this->kMersToGenerate[i].position));
            }
        }

    private:
        const std::vector<KMerPositionInfo> &kMersToGenerate;
        KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig;
        std::vector<KMerStringCreatorForSequence<input_vector_t, input_elem_t>> kmerStringCreators;
        std::function<std::string(int, int)> createKMerFunc;
        std::vector<int> gapsAccumulated;
        std::vector<std::string> &resultStrings;
        int resultOffset;

        inline void prepareKMerStringsCreators() {
            this->kmerStringCreators.reserve(kMerTaskConfig.sequencesNum);
            for (int i = 0; i < kMerTaskConfig.sequencesNum; ++i) {
                auto seq = std::move(kMerTaskConfig.sequenceGetter(i));
                this->kmerStringCreators.emplace_back(
                        std::move(seq),
                        kMerTaskConfig.gaps, gapsAccumulated,
                        kMerTaskConfig.kMerItemSeparator, kMerTaskConfig.kMerSectionSeparator,
                        kMerTaskConfig.inputToStringItemConverter);
            }
        }

        inline void prepareCreateKMerFunc() {
            this->createKMerFunc = kMerTaskConfig.positionalKMers ?
                                   static_cast<std::function<std::string(int, int)>>(
                                           [this](int seqNum, int pos) {
                                               return kmerStringCreators[seqNum].getPositional(pos);
                                           }) :
                                   static_cast<std::function<std::string(int, int)>>(
                                           [this](int seqNum, int pos) {
                                               return kmerStringCreators[seqNum].get(pos);
                                           });
        }

    };

    template<class input_vector_t, class input_elem_t>
    class KMerStringCreatorForSequence {
    public:
        KMerStringCreatorForSequence(input_vector_t &&sequence,
                                     const std::vector<int> &gaps,
                                     const std::vector<int> &gapsAccumulated,
                                     const std::string &itemSeparator,
                                     const std::string &sectionSeparator,
                                     InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) :
                sequence(std::move(sequence)),
                gapsAccumulated(gapsAccumulated),
                itemSeparator(itemSeparator),
                sectionSeparator(sectionSeparator),
                inputToStringItemConverter(inputToStringItemConverter) {
            this->kmerInfoSuffix = prepareKMerInfoSuffix(gaps);
        }

        inline std::string get(int begin) const {
            int totalSize = getTotalSize(begin, itemSeparator.size());
            std::string res;
            res.reserve(totalSize);
            res += inputToStringItemConverter(sequence[begin]);
            for (const int &accGap: this->gapsAccumulated) {
                res += itemSeparator + inputToStringItemConverter(sequence[begin + accGap]);
            }
            return kmerInfoSuffix.size() > 0
                   ? res + sectionSeparator + kmerInfoSuffix
                   : res;
        }

        inline std::string getPositional(int begin) const {
            std::string withoutPositionString = this->get(begin);
            return std::to_string(begin + 1) + sectionSeparator + withoutPositionString;
        }

    private:
        input_vector_t sequence;
        std::string itemSeparator;
        std::string sectionSeparator;
        std::string kmerInfoSuffix;
        const std::vector<int> &gapsAccumulated;
        InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter;

        inline std::size_t getTotalSize(int begin, int separatorLength) const {
            return std::accumulate(
                    std::begin(this->gapsAccumulated),
                    std::end(this->gapsAccumulated),
                    inputToStringItemConverter(sequence[begin]).size(),
                    [this, &begin, &separatorLength](int r, int accGap) {
                        return r + inputToStringItemConverter(sequence[begin + accGap]).size() + separatorLength;
                    }
            );
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

}

#endif
