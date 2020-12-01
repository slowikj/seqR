#ifndef SEQR_KMER_TASK_SOLVER_H
#define SEQR_KMER_TASK_SOLVER_H

#include <Rcpp.h>
#include "kmer_task_config.h"
#include "kmer_manager.h"
#include "kmer_counting_common_algorithm.h"
#include "dictionary/supported_dict_names.h"
#include "dictionary/linear_list_dictionary.h"
#include "gapped_kmer_counter.h"
#include "kmer_counter.h"
#include "dictionary/unordered_map_wrapper.h"
#include <functional>
#include <vector>

const std::string DEFAULT_KMER_ITEM_SEPARATOR = ".";

const std::string DEFAULT_KMER_SECTION_SEPARATOR = "_";

template<class input_vector_t,
        class input_elem_t,
        class alphabet_encoding_t,
        template<typename key, typename value> class kmer_dictionary_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   alphabet_encoding_t &alphabetEncoding,
                   std::function<ComplexHasher()> &complexHasherFactory,
                   KMerCountingResult &kMerCountingResult) {
    KMerCountingProc_t<input_vector_t, input_elem_t, alphabet_encoding_t, kmer_dictionary_t> countingProc = [&complexHasherFactory](
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            alphabet_encoding_t &encoding) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                computeKMersForAllSequences<input_vector_t, kmer_dictionary_t>(
                        kMerTaskConfig.sequencesNum,
                        [&kMerTaskConfig, &complexHasherFactory, &encoding]
                                (input_vector_t &v) -> KMerManager<kmer_dictionary_t> {
                            return countKMers<input_vector_t, alphabet_encoding_t, kmer_dictionary_t>(
                                    kMerTaskConfig.k,
                                    v,
                                    encoding,
                                    kMerTaskConfig.positionalKMers,
                                    kMerTaskConfig.withKMerCounts,
                                    std::move(complexHasherFactory())
                            );
                        },
                        kMerTaskConfig.sequenceGetter,
                        kMerTaskConfig.parallelMode
                ));
    };
    updateKMerCountingResult<input_vector_t, input_elem_t, alphabet_encoding_t, kmer_dictionary_t>(
            kMerTaskConfig, alphabetEncoding, countingProc, kMerCountingResult);
}

template<class input_vector_t,
        class input_elem_t,
        class alphabet_encoding_t,
        template<typename key, typename value> class kmer_dictionary_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   alphabet_encoding_t &alphabetEncoding,
                   std::vector<PolynomialSingleHasherConfig> &hasherConfigs,
                   KMerCountingResult &kMerCountingResult) {
    KMerCountingProc_t<input_vector_t, input_elem_t, alphabet_encoding_t, kmer_dictionary_t> countingProc = [&hasherConfigs](
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            alphabet_encoding_t &encoding) -> std::vector<KMerManager<kmer_dictionary_t>> {
        std::size_t totalKMerSize = getTotalKMerSize(kMerTaskConfig.gaps);
        return std::move(
                computeKMersForAllSequences<input_vector_t, kmer_dictionary_t>(
                        kMerTaskConfig.sequencesNum,
                        [&kMerTaskConfig, &encoding, &totalKMerSize, &hasherConfigs]
                                (input_vector_t &v) -> KMerManager<kmer_dictionary_t> {
                            return countGappedKMers<input_vector_t, alphabet_encoding_t, kmer_dictionary_t>(
                                    kMerTaskConfig.gaps,
                                    totalKMerSize,
                                    v,
                                    encoding,
                                    kMerTaskConfig.positionalKMers,
                                    kMerTaskConfig.withKMerCounts,
                                    hasherConfigs
                            );
                        },
                        kMerTaskConfig.sequenceGetter,
                        kMerTaskConfig.parallelMode
                ));
    };
    updateKMerCountingResult<input_vector_t, input_elem_t, alphabet_encoding_t, kmer_dictionary_t>(
            kMerTaskConfig, alphabetEncoding, countingProc, kMerCountingResult);
}

template<class input_vector_t,
        class input_elem_t,
        class alphabet_encoding_t,
        class algorithm_params_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   alphabet_encoding_t &alphabetEncoding,
                   const std::string &kmerDictionaryName,
                   algorithm_params_t &algorithmParams,
                   KMerCountingResult &kMerCountingResult) {
    if (kmerDictionaryName == UNORDERED_MAP_NAME) {
        computeResult<input_vector_t, input_elem_t, alphabet_encoding_t, UnorderedMapWrapper>(
                kMerTaskConfig, alphabetEncoding, algorithmParams, kMerCountingResult);

    } else if (kmerDictionaryName == LINEAR_LIST_NAME) {
        computeResult<input_vector_t, input_elem_t, alphabet_encoding_t, LinearListDictionary>(
                kMerTaskConfig, alphabetEncoding, algorithmParams, kMerCountingResult);

    } else {
        std::string errorMessage = "unsupported k-mer dictionary name: " + kmerDictionaryName;
        throw Rcpp::exception(errorMessage.c_str());
    }
}

inline void printLogIfVerbose(bool verbose, int begin, int end) {
    if (verbose) {
        Rcpp::Rcout << "Start processing sequences (batch: [" << begin + 1 << "-" << end << "])..." << std::endl;
    }
}

inline
Rcpp::List computeKMersInBatches(const std::function<void(KMerCountingResult &, int, int)> &batchFunc,
                                 int sequencesNum,
                                 int batchSize,
                                 bool verbose) {
    KMerCountingResult kMerCountingResult;
    for (int begin = 0; begin < sequencesNum; begin += batchSize) {
        int end = std::min(begin + batchSize, sequencesNum);
        printLogIfVerbose(verbose, begin, end);
        Rcpp::checkUserInterrupt();
        batchFunc(kMerCountingResult, begin, end);
    }
    return kMerCountingResult.toRcppList();
}

#endif //SEQR_KMER_TASK_SOLVER_H
