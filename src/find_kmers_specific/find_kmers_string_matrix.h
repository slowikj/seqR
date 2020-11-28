//
// Created by slowik on 28.11.2020.
//

#ifndef SEQR_FIND_KMERS_STRING_MATRIX_H
#define SEQR_FIND_KMERS_STRING_MATRIX_H

#include <Rcpp.h>
#include <vector>
#include "../kmer_task_config.h"
#include "../default_alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::StringMatrix &sequenceMatrix,
                             Rcpp::StringVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             bool verbose,
                             algorithm_params_t &algorithmParams) {
    auto cppAlphabet = std::move(Rcpp::as<std::vector<std::string>>(alphabet));
    auto alphabetEncoding = std::move(
            getDefaultAlphabetEncoder<std::vector<std::string>, std::string, short, UnorderedMapWrapper>(cppAlphabet));

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
                DefaultAlphabetEncoder<std::string, short, UnorderedMapWrapper>,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize, verbose);
}


#endif //SEQR_FIND_KMERS_STRING_MATRIX_H
