#ifndef SEQR_FIND_KMERS_NUMERIC_MATRIX_H
#define SEQR_FIND_KMERS_NUMERIC_MATRIX_H

#include <Rcpp.h>
#include <vector>
#include "../kmer_task_config.h"
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
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

template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::NumericMatrix &sequenceMatrix,
                             Rcpp::NumericVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             bool verbose,
                             algorithm_params_t &algorithmParams) {
    using encoded_elem_t = config::encoded_elem_t;
    auto alphabetEncoding = std::move(
            getDefaultAlphabetEncoder<Rcpp::NumericVector, double, encoded_elem_t, UnorderedMapWrapper>(alphabet));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<RcppParallel::RMatrix<double>::Row, double> kMerTaskConfig(
                (seqEnd - seqBegin),
                getRMatrixRowGetter<Rcpp::NumericMatrix, decltype(alphabetEncoding)::input_elem_t>(sequenceMatrix,
                                                                                                   seqBegin),
                gaps,
                positionalKMers,
                withKMerCounts,
                getDoubleToStringConverter(3),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                RcppParallel::RMatrix<double>::Row,
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

#endif //SEQR_FIND_KMERS_NUMERIC_MATRIX_H
