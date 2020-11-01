#ifndef SEQR_KMER_TASK_SOLVER_TYPE_SPECIFIC_H
#define SEQR_KMER_TASK_SOLVER_TYPE_SPECIFIC_H


#include <Rcpp.h>
#include <vector>
#include "kmer_task_config.h"
#include "alphabet_encoder.h"
#include "dictionary/unordered_map_wrapper.h"
#include "kmer_counting_result.h"
#include "kmer_task_solver.h"

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::StringMatrix &sequenceMatrix,
                             Rcpp::StringVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    auto cppAlphabet = std::move(Rcpp::as<std::vector<std::string>>(alphabet));
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<std::vector<std::string>, std::string, short, UnorderedMapWrapper>(cppAlphabet));

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
                             Rcpp::IntegerVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<Rcpp::IntegerVector, int, short, UnorderedMapWrapper>(alphabet));

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
                             Rcpp::NumericVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<Rcpp::NumericVector, double, short, UnorderedMapWrapper>(alphabet));

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


#endif //SEQR_KMER_TASK_SOLVER_TYPE_SPECIFIC_H
