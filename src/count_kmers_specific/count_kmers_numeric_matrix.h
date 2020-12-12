#ifndef SEQR_COUNT_KMERS_NUMERIC_MATRIX_H
#define SEQR_COUNT_KMERS_NUMERIC_MATRIX_H

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/stl_unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"

inline InputToStringItemConverter_t<double> getDoubleToStringConverter(int decimalPrecision) {
    return [decimalPrecision](const double &elem) -> std::string {
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(decimalPrecision);
        stream << elem;
        return stream.str();
    };
}

template<class algorithm_params_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::List commonCountKMersSpecific(Rcpp::NumericMatrix &sequenceMatrix,
                                    Rcpp::NumericVector &alphabet,
                                    const UserParams &userParams,
                                    algorithm_params_t &algorithmParams) {
    using encoded_elem_t = config::encoded_elem_t;
    auto alphabetEncoding = std::move(
            alphabetEncoding::getDefaultAlphabetEncoder<Rcpp::NumericVector, double, encoded_elem_t, dictionary::StlUnorderedMapWrapper>(
                    alphabet));

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<RcppParallel::RMatrix<double>::Row, double> kMerTaskConfig(
                (seqEnd - seqBegin),
                getRMatrixRowGetter<Rcpp::NumericMatrix, decltype(alphabetEncoding)::input_elem_t>(
                        sequenceMatrix, seqBegin),
                getDoubleToStringConverter(3),
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR,
                userParams);
        updateKMerCountingResult<
                RcppParallel::RMatrix<double>::Row,
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
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::List parallelCountKMersSpecific(Rcpp::NumericMatrix &sequenceMatrix,
                                      Rcpp::NumericVector &alphabet,
                                      const UserParams &userParams,
                                      algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequenceMatrix, alphabet, userParams, algorithmParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::List sequentialCountKMersSpecific(Rcpp::NumericMatrix &sequenceMatrix,
                                        Rcpp::NumericVector &alphabet,
                                        const UserParams &userParams,
                                        algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequenceMatrix, alphabet, userParams, algorithmParams);
}

#endif //SEQR_COUNT_KMERS_NUMERIC_MATRIX_H
