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
#include "tidysq_encoded_sequence.h"
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

template<class alphabet_t,
        class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        class algorithm_params_t>
inline
void computeResult(alphabet_t &alphabet,
                   KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                   const std::string &kmerDictionaryName,
                   algorithm_params_t &algorithmParams,
                   KMerCountingResult &kMerCountingResult) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabet));
    computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, algorithm_params_t>(
            kMerTaskConfig,
            alphabetEncoding,
            kmerDictionaryName,
            algorithmParams,
            kMerCountingResult
    );
}

template<class algorithm_params_t>
inline
void findKMersSpecific(Rcpp::StringMatrix &sequenceMatrix,
                       int seqBegin,
                       int seqEnd,
                       std::vector<std::string> &alphabet,
                       std::vector<int> &gaps,
                       bool positionalKMers,
                       bool withKMerCounts,
                       const std::string &kmerDictionaryName,
                       algorithm_params_t &algorithmParams,
                       KMerCountingResult &kMerCountingResult) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    KMerTaskConfig<SafeMatrixSequenceWrapper<std::string>::Row, std::string> kMerTaskConfig(
            (seqEnd - seqBegin),
            getSafeMatrixRowGetter<std::string>(safeMatrixWrapper, seqBegin),
            gaps,
            positionalKMers,
            withKMerCounts,
            getStringToStringConverter(),
            DEFAULT_KMER_ITEM_SEPARATOR,
            DEFAULT_KMER_SECTION_SEPARATOR);
    computeResult<
            std::vector<std::string>,
            SafeMatrixSequenceWrapper<std::string>::Row,
            std::string,
            short,
            UnorderedMapWrapper,
            algorithm_params_t>(alphabet,
                                kMerTaskConfig,
                                kmerDictionaryName,
                                algorithmParams,
                                kMerCountingResult);
}

template<class algorithm_params_t>
inline
void findKMersSpecific(Rcpp::IntegerMatrix &sequenceMatrix,
                       int seqBegin,
                       int seqEnd,
                       std::vector<int> &alphabet,
                       std::vector<int> &gaps,
                       bool positionalKMers,
                       bool withKMerCounts,
                       const std::string &kmerDictionaryName,
                       algorithm_params_t &algorithmParams,
                       KMerCountingResult &kMerCountingResult) {
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
            std::vector<int>,
            RcppParallel::RMatrix<int>::Row,
            int,
            short,
            UnorderedMapWrapper,
            algorithm_params_t>(alphabet,
                                kMerTaskConfig,
                                kmerDictionaryName,
                                algorithmParams,
                                kMerCountingResult);
}

template<class algorithm_params_t>
inline
void findKMersSpecific(Rcpp::NumericMatrix &sequenceMatrix,
                       int seqBegin,
                       int seqEnd,
                       std::vector<double> &alphabet,
                       std::vector<int> &gaps,
                       bool positionalKMers,
                       bool withKMerCounts,
                       const std::string &kmerDictionaryName,
                       algorithm_params_t &algorithmParams,
                       KMerCountingResult &kMerCountingResult) {
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
            std::vector<double>,
            RcppParallel::RMatrix<double>::Row,
            double,
            short,
            UnorderedMapWrapper,
            algorithm_params_t>(alphabet,
                                kMerTaskConfig,
                                kmerDictionaryName,
                                algorithmParams,
                                kMerCountingResult);
}

template<class algorithm_params_t>
inline
void findKMersSpecific(Rcpp::List &sq,
                       int seqBegin,
                       int seqEnd,
                       std::vector<std::string> &alphabet,
                       std::vector<int> &gaps,
                       bool positionalKMers,
                       bool withKMerCounts,
                       const std::string &kmerDictionaryName,
                       algorithm_params_t &algorithmParams,
                       KMerCountingResult &kMerCountingResult) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    std::vector<std::string> safeElementsEncoding = convertRcppVector<std::string, Rcpp::StringVector>(
            elementsEncoding);
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<std::vector<std::string>, unsigned char, UnorderedMapWrapper>(
            alphabet,
            safeElementsEncoding
    ));
    auto encodedSequences = getEncodedTidysqSequences(sq, seqBegin, seqEnd);
    KMerTaskConfig<RcppParallel::RVector<unsigned char>, unsigned char> kMerTaskConfig(
            (seqEnd - seqBegin),
            getTidysqRVectorGetter(encodedSequences),
            gaps,
            positionalKMers,
            withKMerCounts,
            getEncodedTidySqItemToStringConverter(safeElementsEncoding),
            DEFAULT_KMER_ITEM_SEPARATOR,
            DEFAULT_KMER_SECTION_SEPARATOR);
    computeResult<
            RcppParallel::RVector<unsigned char>,
            unsigned char,
            unsigned char,
            UnorderedMapWrapper,
            algorithm_params_t>(kMerTaskConfig,
                                alphabetEncoding,
                                kmerDictionaryName,
                                algorithmParams,
                                kMerCountingResult);
}

template<class seq_t,
        class alphabet_t,
        class algorithm_params_t>
inline
Rcpp::List findKMers(seq_t &sequences,
                     int sequencesNum,
                     alphabet_t &alphabet,
                     std::vector<int> &gaps,
                     bool positionalKMers,
                     bool withKMerCounts,
                     const std::string &kmerDictionaryName,
                     algorithm_params_t &algorithmParams,
                     int batchSize) {
    KMerCountingResult kMerCountingResult;
    for (int seqBegin = 0; seqBegin < sequencesNum; seqBegin += batchSize) {
        int seqEnd = std::min(seqBegin + batchSize, sequencesNum);
        findKMersSpecific<algorithm_params_t>(sequences, seqBegin, seqEnd, alphabet, gaps, positionalKMers,
                                              withKMerCounts,
                                              kmerDictionaryName, algorithmParams, kMerCountingResult);
    }
    return kMerCountingResult.toRcppList();
}

#endif //SEQR_KMER_TASK_SOLVER_H
