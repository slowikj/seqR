#ifndef SEQR_FIND_KMERS_INTEGER_MATRIX_H
#define SEQR_FIND_KMERS_INTEGER_MATRIX_H

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "../common_config.h"

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
                             bool parallelMode,
                             algorithm_params_t &algorithmParams) {
    using encoded_elem_t = config::encoded_elem_t ;
    auto alphabetEncoding = std::move(
            alphabetEncoding::getDefaultAlphabetEncoder<Rcpp::IntegerVector, int, encoded_elem_t, dictionary::UnorderedMapWrapper>(alphabet));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<RcppParallel::RMatrix<int>::Row, int> kMerTaskConfig(
                (seqEnd - seqBegin),
                getRMatrixRowGetter<Rcpp::IntegerMatrix, decltype(alphabetEncoding)::input_elem_t>(sequenceMatrix, seqBegin),
                gaps,
                positionalKMers,
                withKMerCounts,
                parallelMode,
                getIntToStringConverter(),
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR);
        updateKMerCountingResult<
                RcppParallel::RMatrix<int>::Row,
                decltype(alphabetEncoding)::input_elem_t,
                decltype(alphabetEncoding),
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize, verbose);
}

#endif //SEQR_FIND_KMERS_INTEGER_MATRIX_H
