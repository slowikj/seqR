#ifndef SEQR_KMER_TASK_SOLVER_H
#define SEQR_KMER_TASK_SOLVER_H

#include <Rcpp.h>
#include <functional>
#include <vector>
#include "kmer_task_config.h"
#include "alphabet_encoder.h"
#include "kmer_counting_common_algorithm.h"
#include "dictionary/supported_dict_names.h"
#include "dictionary/linear_list_dictionary.h"
#include "gapped_kmer_counter.h"
#include "kmer_counter.h"
#include "kmer_counting_result.h"
#include "dictionary/unordered_map_wrapper.h"

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> &parallelKMerCountingProc,
                   KMerCountingResult &kMerCountingResult) {
    parallelComputeKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
            kMerTaskConfig,
            alphabetEncoding,
            parallelKMerCountingProc,
            kMerCountingResult);
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value> class kmer_dictionary_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   std::function<ComplexHasher()> &complexHasherFactory,
                   KMerCountingResult &kMerCountingResult) {
    ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> countingProc = [&complexHasherFactory](
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                parallelComputeKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                        kMerTaskConfig, encoding, complexHasherFactory));
    };
    computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
            kMerTaskConfig, alphabetEncoding, countingProc, kMerCountingResult);
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value> class kmer_dictionary_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   std::vector<PolynomialSingleHasherConfig> &hasherConfigs,
                   KMerCountingResult &kMerCountingResult) {
    ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> countingProc = [&hasherConfigs](
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                parallelComputeGappedKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                        kMerTaskConfig, encoding, hasherConfigs));
    };
    computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
            kMerTaskConfig, alphabetEncoding, countingProc, kMerCountingResult);
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        class algorithm_params_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   const std::string &kmerDictionaryName,
                   algorithm_params_t &algorithmParams,
                   KMerCountingResult &kMerCountingResult) {
    if (kmerDictionaryName == UNORDERED_MAP_NAME) {
        computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, UnorderedMapWrapper>(
                kMerTaskConfig, alphabetEncoding, algorithmParams, kMerCountingResult);

    } else if (kmerDictionaryName == LINEAR_LIST_NAME) {
        computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, LinearListDictionary>(
                kMerTaskConfig, alphabetEncoding, algorithmParams, kMerCountingResult);
    } else {
        throw std::runtime_error("unsupported k-mer dictionary name " + kmerDictionaryName);
    }
}

#endif //SEQR_KMER_TASK_SOLVER_H
