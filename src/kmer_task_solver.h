#ifndef SEQR_KMER_TASK_SOLVER_H
#define SEQR_KMER_TASK_SOLVER_H

#include <Rcpp.h>
#include "kmer_task_config.h"
#include "alphabet_encoder.h"
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
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
void computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> &parallelKMerCountingProc,
                   KMerCountingResult &kMerCountingResult) {
    parallelGetKMerCounts<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
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
        std::string errorMessage = "unsupported k-mer dictionary name: " + kmerDictionaryName;
        throw Rcpp::exception(errorMessage.c_str());
    }
}

inline
Rcpp::List computeKMersInBatches(const std::function<void(KMerCountingResult &, int, int)> &batchFunc,
                                 int sequencesNum,
                                 int batchSize) {
    KMerCountingResult kMerCountingResult;
    for (int begin = 0; begin < sequencesNum; begin += batchSize) {
        int end = std::min(begin + batchSize, sequencesNum);
        batchFunc(kMerCountingResult, begin, end);
    }
    return kMerCountingResult.toRcppList();
}

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::StringMatrix &sequenceMatrix,
                             std::vector<std::string> &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<std::vector<std::string>, std::string, short, UnorderedMapWrapper>(alphabet));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        SafeSequencesMatrixWrapper<std::string> safeMatrixWrapper(sequenceMatrix, seqBegin, seqEnd);
        KMerTaskConfig<SafeSequencesMatrixWrapper<std::string>::Row, std::string> kMerTaskConfig(
                (seqEnd - seqBegin),
                getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
                gaps,
                positionalKMers,
                withKMerCounts,
                getStringToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                SafeSequencesMatrixWrapper<std::string>::Row,
                std::string,
                short,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize);
}

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::IntegerMatrix &sequenceMatrix,
                             std::vector<int> &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<std::vector<int>, int, short, UnorderedMapWrapper>(alphabet));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<RcppParallel::RMatrix<int>::Row, int> kMerTaskConfig(
                (seqEnd - seqBegin),
                getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix, seqBegin),
                gaps,
                positionalKMers,
                withKMerCounts,
                getIntToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                RcppParallel::RMatrix<int>::Row,
                int,
                short,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize);
}

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::NumericMatrix &sequenceMatrix,
                             std::vector<double> &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<std::vector<double>, double, short, UnorderedMapWrapper>(alphabet));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<RcppParallel::RMatrix<double>::Row, double> kMerTaskConfig(
                (seqEnd - seqBegin),
                getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix, seqBegin),
                gaps,
                positionalKMers,
                withKMerCounts,
                getDoubleToStringConverter(3),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                RcppParallel::RMatrix<double>::Row,
                double,
                short,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize);
}

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::List &sequences,
                             Rcpp::StringVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    std::string alphabetStr("", alphabet.size());
    for (int i = 0; i < alphabet.size(); ++i) {
        alphabetStr[i] = Rcpp::as<char>(alphabet[i]);
    }
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<std::string, char, short, UnorderedMapWrapper>(alphabetStr));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        SafeSequencesStringListWrapper sequenceWrapper(sequences, seqBegin, seqEnd);
        KMerTaskConfig<SafeSequencesStringListWrapper::Row, char> kMerTaskConfig(
                (seqEnd - seqBegin),
                getStringSequenceGetter(sequenceWrapper),
                gaps,
                positionalKMers,
                withKMerCounts,
                getCharToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                SafeSequencesStringListWrapper::Row,
                char,
                short,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequences.size(), batchSize);
}

#endif //SEQR_KMER_TASK_SOLVER_H
