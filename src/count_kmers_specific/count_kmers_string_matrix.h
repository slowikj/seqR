#pragma once

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/stl_unordered_map_wrapper.h"
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
                        encoded_elem_t notAllowedEncodingNum,
                        int begin, int end) :
            nrow(end - begin),
            ncol(matrix.ncol()),
            encoded(nrow * ncol) {
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
        return Row(encoded, index * ncol, ncol);
    }

private:
    int ncol, nrow;
    std::vector<encoded_elem_t> encoded;
};

template<class algorithm_params_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::List parallelCountKMersSpecific(Rcpp::StringMatrix &sequenceMatrix,
                                      Rcpp::StringVector &alphabet,
                                      const UserParams &userParams,
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

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult, int seqBegin, int seqEnd) {
        EncodedStringMatrix<encoded_elem_t> safeMatrixWrapper(sequenceMatrix,
                                                              stringAlphabetEncoder,
                                                              alphabetBeginCnt - 1,
                                                              seqBegin, seqEnd);
        KMerTaskConfig<typename decltype(safeMatrixWrapper)::Row, decltype(alphabetEncoder)::input_elem_t> kMerTaskConfig(
                (seqEnd - seqBegin),
                [&safeMatrixWrapper](int rowNum) -> typename EncodedStringMatrix<encoded_elem_t>::Row {
                    return safeMatrixWrapper.row(rowNum);
                },
                [&alphabetStrings](const encoded_elem_t &encodedElem) -> std::string {
                    return alphabetStrings[encodedElem - 1];
                },
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR,
                userParams);
        updateKMerCountingResult<
                typename decltype(safeMatrixWrapper)::Row,
                decltype(alphabetEncoder)::input_elem_t,
                decltype(alphabetEncoder),
                kmer_dictionary_t>(kMerTaskConfig,
                                   alphabetEncoder,
                                   algorithmParams,
                                   kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sequenceMatrix.nrow(), userParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::List sequentialCountKMersSpecific(Rcpp::StringMatrix &sequenceMatrix,
                                        Rcpp::StringVector &alphabet,
                                        const UserParams &userParams,
                                        algorithm_params_t &algorithmParams) {
    using encoded_elem_t = config::encoded_elem_t;
    auto alphabetEncoder = alphabetEncoding::getDefaultAlphabetEncoder<
            Rcpp::StringVector, Rcpp::StringVector::stored_type, encoded_elem_t, dictionary::StlUnorderedMapWrapper>(
            alphabet);

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult, int seqBegin, int seqEnd) {
        KMerTaskConfig<typename Rcpp::StringMatrix::Row, decltype(alphabetEncoder)::input_elem_t> kMerTaskConfig(
                (seqEnd - seqBegin),
                [&](int index) -> typename Rcpp::StringMatrix::Row {
                    return sequenceMatrix(index + seqBegin, Rcpp::_);
                },
                [](const Rcpp::StringMatrix::Row::elem_type &elem) -> std::string {
                    return Rcpp::as<std::string>(elem);
                },
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR,
                userParams);
        updateKMerCountingResult<
                typename Rcpp::StringMatrix::Row,
                decltype(alphabetEncoder)::input_elem_t,
                decltype(alphabetEncoder),
                kmer_dictionary_t>(kMerTaskConfig,
                                   alphabetEncoder,
                                   algorithmParams,
                                   kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sequenceMatrix.nrow(), userParams);
}

