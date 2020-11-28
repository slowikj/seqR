//
// Created by slowik on 28.11.2020.
//

#ifndef SEQR_FIND_KMERS_INTEGER_MATRIX_H
#define SEQR_FIND_KMERS_INTEGER_MATRIX_H

#include <Rcpp.h>
#include <vector>
#include "../kmer_task_config.h"
#include "../default_alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"

inline InputToStringItemConverter_t<int> getIntToStringConverter() {
    return [](const int &elem) -> std::string {
        return std::to_string(elem);
    };
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
                             bool verbose,
                             algorithm_params_t &algorithmParams) {
    auto alphabetEncoding = std::move(
            getDefaultAlphabetEncoder<Rcpp::IntegerVector, int, short, UnorderedMapWrapper>(alphabet));

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
                DefaultAlphabetEncoder<int, short, UnorderedMapWrapper>,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize, verbose);
}

#endif //SEQR_FIND_KMERS_INTEGER_MATRIX_H
