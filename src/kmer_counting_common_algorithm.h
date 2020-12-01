#ifndef KMER_COUNTING_COMMON_H
#define KMER_COUNTING_COMMON_H

#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include <vector>
#include <memory>
#include <functional>
#include "sequence_getter.h"
#include "kmer_manager.h"
#include "kmer_strings_creator.h"
#include "input_to_string_item_converter.h"
#include "kmer_counting_result.h"
#include "kmer_position_info.h"
#include "kmer_task_config.h"

template<class input_vector_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
using CountingKMersProc_t = std::function<KMerManager<kmer_dictionary_t>(input_vector_t &)>;

template<class input_vector_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
class KMerCounterWorker : public RcppParallel::Worker {
public:
    KMerCounterWorker(int rowsNum,
                      CountingKMersProc_t<input_vector_t, kmer_dictionary_t> countingKMersProc,
                      SequenceGetter_t<input_vector_t> sequenceGetter) :
            countingKMersProc(countingKMersProc),
            sequenceGetter(sequenceGetter) {
        kMers.resize(rowsNum);
    }

    inline void operator()(size_t begin, size_t end) override {
        for (int rowNum = begin; rowNum < end; ++rowNum) {
            auto row = std::move(sequenceGetter(rowNum));
            kMers[rowNum] = std::move(countingKMersProc(row));
        }
    }

private:
    CountingKMersProc_t<input_vector_t, kmer_dictionary_t> countingKMersProc;
    SequenceGetter_t<input_vector_t> sequenceGetter;

public:
    std::vector<KMerManager<kmer_dictionary_t>> kMers;
};

template<class input_vector_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
std::vector<KMerManager<kmer_dictionary_t>> computeKMersForAllSequences(
        int rowsNum,
        CountingKMersProc_t<input_vector_t, kmer_dictionary_t> countingProc,
        SequenceGetter_t<input_vector_t> sequenceGetter,
        bool parallelMode) {
    KMerCounterWorker<input_vector_t, kmer_dictionary_t> worker(rowsNum, countingProc, sequenceGetter);
    if (parallelMode) {
        RcppParallel::parallelFor(0, rowsNum, worker);
    } else {
        worker(0, rowsNum);
    }
    return std::move(worker.kMers);
}

template<class input_vector_t, class input_elem_t,
        class alphabet_encoding_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
using KMerCountingProc_t = std::function<std::vector<KMerManager<kmer_dictionary_t>>(
        KMerTaskConfig<input_vector_t, input_elem_t> &,
        alphabet_encoding_t &)>;

template<class input_vector_t, class input_elem_t,
        class alphabet_encoding_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
void updateKMerCountingResult(
        KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConf,
        alphabet_encoding_t &alphabetEncoding,
        KMerCountingProc_t<input_vector_t, input_elem_t, alphabet_encoding_t, kmer_dictionary_t> kMerCountingProc,
        KMerCountingResult &kMerCountingResult) {
    auto kMersManagers = std::move(kMerCountingProc(kMerTaskConf, alphabetEncoding));

    std::vector<KMerPositionInfo> kMersToCreate;
    for (int seqNum = 0; seqNum < kMersManagers.size(); ++seqNum) {
        for (const auto &kMerPair: kMersManagers[seqNum].getDictionary()) {
            bool kMerStringNeedsCreation = kMerCountingResult.addKMer(kMerPair.first, seqNum, kMerPair.second.cnt);
            if (kMerStringNeedsCreation) {
                kMersToCreate.emplace_back(seqNum, kMerPair.second.seqStartPosition);
            }
        }
    }

    generateKMerStrings<input_vector_t, input_elem_t>(kMersToCreate,
                                                      kMerTaskConf,
                                                      kMerCountingResult.kMerStrings);

    kMerCountingResult.increaseProcessSequencesNum(kMerTaskConf.sequencesNum);
}

#endif
