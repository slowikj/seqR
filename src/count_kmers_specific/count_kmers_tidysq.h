#ifndef SEQR_COUNT_KMERS_TIDYSQ_H
#define SEQR_COUNT_KMERS_TIDYSQ_H

#include <Rcpp/vector/instantiation.h>
#include "../user_params.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include "../kmer_counting_result.h"
#include "../kmer_task_config.h"
#include "../kmer_task_solver.h"
#include "../common_config.h"
#include <vector>

class TidysqAlphabetEncoder {
public:
    using input_elem_t = tidysq::LetterValue;
    using encoded_elem_t = uint16_t;

    explicit TidysqAlphabetEncoder(
            const tidysq::Alphabet &sqAlphabet,
            const Rcpp::StringVector &kMerAlphabet) :
            size_(kMerAlphabet.size()) {
        prepareIsAllowedArray(sqAlphabet, kMerAlphabet);
    }

    [[nodiscard]] inline encoded_elem_t encode(const input_elem_t &inputElem) const {
        return static_cast<encoded_elem_t>(inputElem);
    }

    [[nodiscard]] inline encoded_elem_t encodeUnsafe(const input_elem_t &inputElem) const {
        return encode(inputElem);
    }

    [[nodiscard]] inline bool isAllowed(const input_elem_t &inputElem) const {
        return isAllowed_[inputElem];
    }

    [[nodiscard]] inline std::size_t size() const {
        return size_;
    }

private:
    std::array<bool, USHRT_MAX> isAllowed_{};
    std::size_t size_;

    void prepareIsAllowedArray(
            const tidysq::Alphabet &sqAlphabet,
            const Rcpp::StringVector &kMerAlphabet) {
        for (const auto &kMerAlphabetElem: kMerAlphabet) {
            auto elem = Rcpp::as<std::string>(kMerAlphabetElem);
            isAllowed_[sqAlphabet.match_value(elem)] = true;
        }
    }
};

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List commonCountKMersSpecific(Rcpp::List &sequences,
                                    Rcpp::StringVector &kMerAlphabet,
                                    const UserParams &userParams,
                                    algorithm_params_t &algorithmParams) {
    auto sq = tidysq::import_from_R(sequences, tidysq::constants::DEFAULT_NA_LETTER);
    auto sqAlphabet = sq.alphabet();
    TidysqAlphabetEncoder alphabetEncoder(sqAlphabet, kMerAlphabet);

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult, int seqBegin, int seqEnd) {
        auto sqUnpackedInts = sq.unpack<tidysq::STD_IT, tidysq::INTS_PT>(seqBegin, seqEnd);

        KMerTaskConfig<decltype(sqUnpackedInts)::ElementType, decltype(alphabetEncoder)::input_elem_t> kMerTaskConfig(
                (seqEnd - seqBegin),
                [&sqUnpackedInts, seqBegin](int index) -> decltype(sqUnpackedInts)::ElementType {
                    return sqUnpackedInts[seqBegin + index];
                },
                [&sqAlphabet](const tidysq::LetterValue &elem) -> std::string { return sqAlphabet[elem]; },
                config::DEFAULT_KMER_ITEM_SEPARATOR,
                config::DEFAULT_KMER_SECTION_SEPARATOR,
                userParams);
        updateKMerCountingResult<
                decltype(sqUnpackedInts)::ElementType,
                decltype(alphabetEncoder)::input_elem_t,
                decltype(alphabetEncoder),
                kmer_dictionary_t>(kMerTaskConfig,
                                   alphabetEncoder,
                                   algorithmParams,
                                   kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sq.size(), userParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List parallelCountKMersSpecific(Rcpp::List &sequences,
                                      Rcpp::StringVector &alphabet,
                                      const UserParams &userParams,
                                      algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequences, alphabet, userParams, algorithmParams);
}

template<class algorithm_params_t,
        template<typename key, typename value, class...> class kmer_dictionary_t>
inline
Rcpp::List sequentialCountKMersSpecific(Rcpp::List &sequences,
                                        Rcpp::StringVector &alphabet,
                                        const UserParams &userParams,
                                        algorithm_params_t &algorithmParams) {
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
            sequences, alphabet, userParams, algorithmParams);
}

#endif //SEQR_COUNT_KMERS_TIDYSQ_H
