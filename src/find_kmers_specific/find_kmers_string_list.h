#ifndef SEQR_FIND_KMERS_STRING_LIST_H
#define SEQR_FIND_KMERS_STRING_LIST_H

#include <Rcpp.h>
#include <vector>
#include "../alphabet_encoder/default_alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "safe_sequences_wrapper.h"
#include <limits>

class SafeStringListWrapper : public BaseSequencesWrapper<std::string, char> {
public:
    explicit SafeStringListWrapper(const Rcpp::List &inputVector)
            : SafeStringListWrapper(inputVector, 0, inputVector.size()) {}

    SafeStringListWrapper(const Rcpp::List &inputVector, std::size_t begin, std::size_t end) {
        this->initSequences(inputVector, begin, end);
    }

private:
    void initSequences(const Rcpp::List &inputList, int begin, int end) {
        this->sequences_.resize(end - begin);
        for (int i = begin; i < end; ++i) {
            Rcpp::StringVector listElem = inputList[i];
            this->sequences_[i - begin] = std::move(Rcpp::as<std::string>(listElem[0]));
        }
    }
};

class StringListAlphabetEncoder {
public:
    using encoded_elem_t = char;
    using input_elem_t = char;

    explicit StringListAlphabetEncoder(std::array<bool, CHAR_MAX> &isElementAllowed, std::size_t size) :
            isElementAllowed(isElementAllowed), size_(size) {
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
    std::array<bool, CHAR_MAX> &isElementAllowed;
    std::size_t size_;
};

inline SequenceGetter_t<typename SafeStringListWrapper::Row>
getStringSequenceGetter(SafeStringListWrapper &sequencesWrapper, int rowOffset = 0) {
    return [&sequencesWrapper, rowOffset](int rowNum) -> typename SafeStringListWrapper::Row {
        return std::move(sequencesWrapper.row(rowNum + rowOffset));
    };
}

inline InputToStringItemConverter_t<char> getCharToStringConverter() {
    return [](const char &elem) -> std::string {
        return std::string(1, elem);
    };
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
                             bool verbose,
                             bool parallelMode,
                             algorithm_params_t &algorithmParams) {
    std::string alphabetStr("", alphabet.size());
    for (int i = 0; i < alphabet.size(); ++i) {
        alphabetStr[i] = Rcpp::as<char>(alphabet[i]);
    }
    std::array<bool, CHAR_MAX> isElementAllowed{};
    for (const char &allowedElem: alphabetStr) {
        isElementAllowed[allowedElem] = true;
    }
    StringListAlphabetEncoder alphabetEncoder(isElementAllowed, alphabet.size());

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        SafeStringListWrapper sequenceWrapper(sequences, seqBegin, seqEnd);
        KMerTaskConfig<SafeStringListWrapper::Row, decltype(alphabetEncoder)::input_elem_t> kMerTaskConfig(
                (seqEnd - seqBegin),
                getStringSequenceGetter(sequenceWrapper),
                gaps,
                positionalKMers,
                withKMerCounts,
                parallelMode,
                getCharToStringConverter(),
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR);
        updateKMerCountingResult<
                SafeStringListWrapper::Row,
                decltype(alphabetEncoder)::input_elem_t,
                decltype(alphabetEncoder),
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoder,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequences.size(), batchSize, verbose);
}


#endif //SEQR_FIND_KMERS_STRING_LIST_H
