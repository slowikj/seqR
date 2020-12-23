#pragma once

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/stl_unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "../common_config.h"
#include "rmatrix_row_getter.h"

inline InputToStringItemConverter_t<int> getIntToStringConverter() {
    return [](const int &elem) -> std::string {
        return std::to_string(elem);
    };
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List commonCountKMersSpecific(Rcpp::IntegerMatrix &sequenceMatrix,
                                    Rcpp::IntegerVector &alphabet,
                                    const UserParams &userParams,
                                    algorithm_params_t &algorithmParams) {
    using encoded_elem_t = config::encoded_elem_t;
    auto alphabetEncoding = alphabetEncoding::getDefaultAlphabetEncoder<Rcpp::IntegerVector, int, encoded_elem_t, dictionary::StlUnorderedMapWrapper>(
            alphabet);

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<RcppParallel::RMatrix<int>::Row, int> kMerTaskConfig(
                (seqEnd - seqBegin),
                getRMatrixRowGetter<Rcpp::IntegerMatrix, decltype(alphabetEncoding)::input_elem_t>(
                        sequenceMatrix, seqBegin),
                getIntToStringConverter(),
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR,
                userParams);
        updateKMerCountingResult<
                RcppParallel::RMatrix<int>::Row,
                decltype(alphabetEncoding)::input_elem_t,
                decltype(alphabetEncoding),
                kmer_dictionary_t>(kMerTaskConfig,
                                   alphabetEncoding,
                                   algorithmParams,
                                   kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sequenceMatrix.nrow(), userParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List parallelCountKMersSpecific(Rcpp::IntegerMatrix &sequenceMatrix,
                                      Rcpp::IntegerVector &alphabet,
                                      const UserParams &userParams,
                                      algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequenceMatrix, alphabet, userParams, algorithmParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List sequentialCountKMersSpecific(Rcpp::IntegerMatrix &sequenceMatrix,
                                        Rcpp::IntegerVector &alphabet,
                                        const UserParams &userParams,
                                        algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequenceMatrix, alphabet, userParams, algorithmParams);
}

