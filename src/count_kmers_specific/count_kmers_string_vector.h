#pragma once

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/stl_unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "safe_sequences_wrapper.h"
#include <limits>

class SafeStringVectorWrapper : public BaseSequencesWrapper<std::string, char> {
public:
    explicit SafeStringVectorWrapper(const Rcpp::StringVector &inputVector)
            : SafeStringVectorWrapper(inputVector, 0, inputVector.size()) {}

    SafeStringVectorWrapper(const Rcpp::StringVector &inputVector, std::size_t begin, std::size_t end) {
        this->initSequences(inputVector, begin, end);
    }

private:
    void initSequences(const Rcpp::StringVector &inputVector, int begin, int end) {
        this->sequences_.resize(end - begin);
        for (int i = begin; i < end; ++i) {
            this->sequences_[i - begin] = Rcpp::as<std::string>(inputVector[i]);
        }
    }
};

class StringListAlphabetEncoder {
public:
    using encoded_elem_t = char;
    using input_elem_t = char;

    explicit StringListAlphabetEncoder(Rcpp::StringVector &alphabet) :
            size_(alphabet.size()) {
        prepareIsElementAllowed(alphabet);
    }

    inline encoded_elem_t encode(const input_elem_t &inputElem) const {
        return inputElem;
    }

    inline encoded_elem_t encodeUnsafe(const input_elem_t &inputElem) const {
        return encode(inputElem);
    }

    inline bool isAllowed(const input_elem_t &inputElem) const {
        return isElementAllowed[inputElem];
    }

    inline std::size_t size() const {
        return size_;
    }

private:
    std::array<bool, CHAR_MAX> isElementAllowed{};
    std::size_t size_;

    inline void prepareIsElementAllowed(Rcpp::StringVector &alphabet) {
        for (const auto &allowedElem: alphabet) {
            isElementAllowed[Rcpp::as<char>(allowedElem)] = true;
        }
    }
};

inline SequenceGetter_t<typename SafeStringVectorWrapper::Row>
getCppStringSequenceGetter(SafeStringVectorWrapper &sequencesWrapper, int rowOffset = 0) {
    return [&sequencesWrapper, rowOffset](int rowNum) -> typename SafeStringVectorWrapper::Row {
        return sequencesWrapper.row(rowNum + rowOffset);
    };
}

inline InputToStringItemConverter_t<char> getCharToStringConverter() {
    return [](const char &elem) -> std::string {
        return std::string(1, elem);
    };
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List commonCountKMersSpecific(Rcpp::StringVector &sequences,
                                    Rcpp::StringVector &alphabet,
                                    const UserParams &userParams,
                                    algorithm_params_t &algorithmParams) {
    StringListAlphabetEncoder alphabetEncoder(alphabet);

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult, int seqBegin, int seqEnd) {
        SafeStringVectorWrapper sequenceWrapper(sequences, seqBegin, seqEnd);
        KMerTaskConfig<SafeStringVectorWrapper::Row, decltype(alphabetEncoder)::input_elem_t> kMerTaskConfig(
                (seqEnd - seqBegin),
                getCppStringSequenceGetter(sequenceWrapper),
                getCharToStringConverter(),
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR,
                userParams);
        updateKMerCountingResult<
                SafeStringVectorWrapper::Row,
                decltype(alphabetEncoder)::input_elem_t,
                decltype(alphabetEncoder),
                kmer_dictionary_t>(kMerTaskConfig,
                                   alphabetEncoder,
                                   algorithmParams,
                                   kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sequences.size(), userParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List parallelCountKMersSpecific(Rcpp::StringVector &sequences,
                                      Rcpp::StringVector &alphabet,
                                      const UserParams &userParams,
                                      algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequences, alphabet, userParams, algorithmParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List sequentialCountKMersSpecific(Rcpp::StringVector &sequences,
                                        Rcpp::StringVector &alphabet,
                                        const UserParams &userParams,
                                        algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequences, alphabet, userParams, algorithmParams);
}
