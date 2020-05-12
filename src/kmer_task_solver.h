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
Rcpp::IntegerMatrix computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                                  AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                  ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> &parallelKMerCountingProc) {
    return std::move(
            getKMerCountsMatrix<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    kMerTaskConfig,
                    alphabetEncoding,
                    parallelKMerCountingProc));
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value> class kmer_dictionary_t>
inline
Rcpp::IntegerMatrix computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                                  AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                  std::function<ComplexHasher()> &complexHasherFactory) {
    ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> countingProc = [&complexHasherFactory](
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                parallelComputeKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                        kMerTaskConfig, encoding, complexHasherFactory));
    };
    return std::move(
            computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    kMerTaskConfig, alphabetEncoding, countingProc));
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value> class kmer_dictionary_t>
inline
Rcpp::IntegerMatrix computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                                  AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                  std::vector<PolynomialSingleHasherConfig> &hasherConfigs) {
    ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t> countingProc = [&hasherConfigs](
            KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                parallelComputeGappedKMersCounts<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                        kMerTaskConfig, encoding, hasherConfigs));
    };
    return std::move(
            computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    kMerTaskConfig, alphabetEncoding, countingProc));
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        class algorithm_params_t>
inline
Rcpp::IntegerMatrix
computeResult(KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
              AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
              const std::string &kmerDictionaryName,
              algorithm_params_t &algorithmParams) {
    if (kmerDictionaryName == UNORDERED_MAP_NAME) {
        return std::move(
                computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, UnorderedMapWrapper>(
                        kMerTaskConfig, alphabetEncoding, algorithmParams));

    } else if (kmerDictionaryName == LINEAR_LIST_NAME) {
        return std::move(
                computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, LinearListDictionary>(
                        kMerTaskConfig, alphabetEncoding, algorithmParams));

    } else {
        throw std::runtime_error("unsupported k-mer dictionary name " + kmerDictionaryName);
    }
}

template<class alphabet_t,
        class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        class algorithm_params_t>
inline
Rcpp::IntegerMatrix computeResult(alphabet_t &alphabet,
                                  KMerTaskConfig<input_vector_t, input_elem_t> &kMerTaskConfig,
                                  const std::string &kmerDictionaryName,
                                  algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabet
            ));
    return std::move(
            computeResult<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, algorithm_params_t>(
                    kMerTaskConfig,
                    alphabetEncoding,
                    kmerDictionaryName,
                    algorithmParams
            ));
}

template<class algorithm_params_t>
inline
Rcpp::IntegerMatrix findKMers(Rcpp::StringMatrix &sequenceMatrix,
                              Rcpp::StringVector &alphabet,
                              std::vector<int> &gaps,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName,
                              algorithm_params_t &algorithmParams) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    std::vector<std::string> convertedAlphabet = std::move(
            convertRcppVector<std::string, Rcpp::StringVector>(alphabet));
    KMerTaskConfig<SafeMatrixSequenceWrapper<std::string>::Row, std::string> kMerTaskConfig(
            sequenceMatrix.nrow(),
            getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
            gaps,
            positionalKMers,
            withKMerCounts,
            getStringToStringConverter(),
            DEFAULT_KMER_ITEM_SEPARATOR,
            DEFAULT_KMER_SECTION_SEPARATOR);
    return std::move(
            computeResult<
                    std::vector<std::string>,
                    SafeMatrixSequenceWrapper<std::string>::Row,
                    std::string,
                    short,
                    UnorderedMapWrapper,
                    algorithm_params_t>(convertedAlphabet,
                                        kMerTaskConfig,
                                        kmerDictionaryName,
                                        algorithmParams));
}

template<class algorithm_params_t>
inline
Rcpp::IntegerMatrix findKMers(Rcpp::IntegerMatrix &sequenceMatrix,
                              Rcpp::IntegerVector &alphabet,
                              std::vector<int> &gaps,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName,
                              algorithm_params_t &algorithmParams) {
    std::vector<int> convertedAlphabet = std::move(convertRcppVector<int, Rcpp::IntegerVector>(alphabet));
    KMerTaskConfig<RcppParallel::RMatrix<int>::Row, int> kMerTaskConfig(
            sequenceMatrix.nrow(),
            getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix),
            gaps,
            positionalKMers,
            withKMerCounts,
            getIntToStringConverter(),
            DEFAULT_KMER_ITEM_SEPARATOR,
            DEFAULT_KMER_SECTION_SEPARATOR);
    return computeResult<
            std::vector<int>,
            RcppParallel::RMatrix<int>::Row,
            int,
            short,
            UnorderedMapWrapper,
            algorithm_params_t>(convertedAlphabet,
                                kMerTaskConfig,
                                kmerDictionaryName,
                                algorithmParams);
}

template<class algorithm_params_t>
inline
Rcpp::IntegerMatrix findKMers(Rcpp::NumericMatrix &sequenceMatrix,
                              Rcpp::NumericVector &alphabet,
                              std::vector<int> &gaps,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName,
                              algorithm_params_t &algorithmParams) {
    std::vector<double> convertedAlphabet = std::move(convertRcppVector<double, Rcpp::NumericVector>(alphabet));
    KMerTaskConfig<RcppParallel::RMatrix<double>::Row, double> kMerTaskConfig(
            sequenceMatrix.nrow(),
            getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix),
            gaps,
            positionalKMers,
            withKMerCounts,
            getDoubleToStringConverter(3),
            DEFAULT_KMER_ITEM_SEPARATOR,
            DEFAULT_KMER_SECTION_SEPARATOR);
    return computeResult<
            std::vector<double>,
            RcppParallel::RMatrix<double>::Row,
            double,
            short,
            UnorderedMapWrapper,
            algorithm_params_t>(convertedAlphabet,
                                kMerTaskConfig,
                                kmerDictionaryName,
                                algorithmParams);
}

template<class algorithm_params_t>
inline
Rcpp::IntegerMatrix findKMers(Rcpp::List &sq,
                              Rcpp::StringVector &alphabet,
                              std::vector<int> &gaps,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName,
                              algorithm_params_t &algorithmParams) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, UnorderedMapWrapper>(
            alphabet,
            elementsEncoding
    ));
    std::vector<std::string> safeElementsEncoding = convertRcppVector<std::string, Rcpp::StringVector>(
            elementsEncoding);
    auto encodedSequences = getEncodedTidysqSequences(sq);
    KMerTaskConfig<RcppParallel::RVector<unsigned char>, unsigned char> kMerTaskConfig(
            sq.size(),
            getTidysqRVectorGetter(encodedSequences),
            gaps,
            positionalKMers,
            withKMerCounts,
            getEncodedTidySqItemToStringConverter(safeElementsEncoding),
            DEFAULT_KMER_ITEM_SEPARATOR,
            DEFAULT_KMER_SECTION_SEPARATOR);
    return computeResult<
            RcppParallel::RVector<unsigned char>,
            unsigned char,
            unsigned char,
            UnorderedMapWrapper,
            algorithm_params_t>(kMerTaskConfig,
                                alphabetEncoding,
                                kmerDictionaryName,
                                algorithmParams);
}

#endif //SEQR_KMER_TASK_SOLVER_H
