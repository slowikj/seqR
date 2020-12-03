#ifndef SEQR_FIND_KMERS_STRING_MATRIX_H
#define SEQR_FIND_KMERS_STRING_MATRIX_H

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "rmatrix_row_getter.h"
#include "../common_config.h"
#include "../alphabet_encoder/identity_alphabet_encoder.h"
#include "encoded_sequence_row.h"

template<class encoded_elem_t>
class EncodedStringMatrix {
public:
    using Row = EncodedSequenceRow<encoded_elem_t>;

    EncodedStringMatrix(Rcpp::StringMatrix &matrix,
                        std::unordered_map<Rcpp::StringMatrix::stored_type, encoded_elem_t> &encoder,
                        int notAllowedEncodingNum,
                        int begin, int end) :
            nrow(end - begin), ncol(matrix.ncol()) {
        this->encoded = std::move(std::shared_ptr<encoded_elem_t[]>(new encoded_elem_t[nrow * ncol]));

        for (int r = begin; r < end; ++r) {
            for (int c = 0; c < ncol; ++c) {
                int encodedIndex = (r - begin) * ncol + c;
                auto matrixElem = matrix(r, c);
                encoded[encodedIndex] = encoder.find(matrixElem) != encoder.end() ?
                                        encoder[matrixElem] :
                                        notAllowedEncodingNum;
            }
        }
    }

    inline Row row(int index) {
        return std::move(Row(encoded, index * ncol, ncol));
    }

private:
    std::shared_ptr<encoded_elem_t[]> encoded;
    int ncol, nrow;
};

template<class encoded_elem_t>
inline SequenceGetter_t<typename EncodedStringMatrix<encoded_elem_t>::Row>
getFastEncodedMatrixGetter(EncodedStringMatrix<encoded_elem_t> &wrapper, int rowOffset = 0) {
    return [&wrapper, rowOffset](int rowNum) -> typename EncodedStringMatrix<encoded_elem_t>::Row {
        return std::move(wrapper.row(rowNum + rowOffset));
    };
}

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
                             bool parallelMode,
                             algorithm_params_t &algorithmParams) {
    using encoded_elem_t = config::encoded_elem_t;
    std::unordered_map<Rcpp::StringVector::stored_type, encoded_elem_t> stringAlphabetEncoder;
    std::vector<std::string> alphabetStrings;
    encoded_elem_t alphabetBeginCnt = 1;
    encoded_elem_t alphabetCnt = alphabetBeginCnt;
    for (const auto &alphabetElem: alphabet) {
        if (stringAlphabetEncoder.find(alphabetElem) == stringAlphabetEncoder.end()) {
            stringAlphabetEncoder[alphabetElem] = alphabetCnt++;
            alphabetStrings.push_back(Rcpp::as<std::string>(alphabetElem));
        }
    }
    alphabetEncoding::IdentityAlphabetEncoder<encoded_elem_t> alphabetEncoder(alphabetBeginCnt, alphabetCnt - 1);

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        EncodedStringMatrix<encoded_elem_t> safeMatrixWrapper(sequenceMatrix,
                                                              stringAlphabetEncoder,
                                                                   alphabetBeginCnt - 1,
                                                              seqBegin, seqEnd);
        KMerTaskConfig<typename decltype(safeMatrixWrapper)::Row, decltype(alphabetEncoder)::input_elem_t> kMerTaskConfig(
                (seqEnd - seqBegin),
                getFastEncodedMatrixGetter<encoded_elem_t>(safeMatrixWrapper),
                gaps,
                positionalKMers,
                withKMerCounts,
                parallelMode,
                [&alphabetStrings](const encoded_elem_t& encodedElem) -> std::string { return alphabetStrings[encodedElem - 1]; },
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                typename decltype(safeMatrixWrapper)::Row,
                decltype(alphabetEncoder)::input_elem_t,
                decltype(alphabetEncoder),
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoder,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequenceMatrix.nrow(), batchSize, verbose);
}


#endif //SEQR_FIND_KMERS_STRING_MATRIX_H
