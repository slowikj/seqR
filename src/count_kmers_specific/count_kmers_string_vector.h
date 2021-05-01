#pragma once

#include <Rcpp.h>
#include <vector>
#include <limits>
#include <string>
#include <array>
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "../encoded_sequence/encoded_sequence_proxy.h"
#include "../common_config.h"

class EncodedStringVector
{
public:
    using init_elem_t = std::string;
    using encoded_elem_t = char;
    using Entry = EncodedSequenceProxy<EncodedStringVector>;

    EncodedStringVector(
        const std::array<bool, CHAR_MAX> &isAllowed,
        Rcpp::StringVector sequences,
        std::size_t seqBegin,
        std::size_t seqEnd)
        : _isAllowed(isAllowed)
    {
        _encode(sequences, seqBegin, seqEnd);
    }

    EncodedStringVector() = delete;

    inline Entry operator[](std::size_t sequenceNum) const
    {
        return Entry(sequenceNum, *this);
    }

    inline encoded_elem_t getElem(std::size_t sequenceNum,
                                  std::size_t offset) const
    {
        return _encodedSequences[sequenceNum][offset];
    }

    inline init_elem_t decode(std::size_t sequenceNum,
                              std::size_t offset) const
    {
        return std::string(1, _encodedSequences[sequenceNum][offset]);
    }

    inline bool isAllowed(std::size_t sequenceNum,
                          std::size_t offset) const
    {
        return _isAllowed[_encodedSequences[sequenceNum][offset]];
    }

    inline std::size_t getSequenceSize(std::size_t sequenceNum) const
    {
        return _encodedSequences[sequenceNum].size();
    }

    inline std::size_t size() const
    {
        return _encodedSequences.size();
    }

private:
    const std::array<bool, CHAR_MAX> &_isAllowed;
    std::vector<std::string> _encodedSequences;

    void _encode(Rcpp::StringVector sequences, std::size_t seqBegin, std::size_t seqEnd)
    {
        _encodedSequences.resize(seqEnd - seqBegin);
        for (std::size_t i = seqBegin; i < seqEnd; ++i)
        {
            _encodedSequences[i - seqBegin] = Rcpp::as<std::string>(sequences[i]);
        }
    }
};

template <class algorithm_params_t,
          template <typename key, typename value, class...> class kmer_dictionary_t>
inline Rcpp::List commonCountKMersSpecific(Rcpp::StringVector &sequences,
                                           Rcpp::StringVector &alphabet,
                                           const UserParams &userParams,
                                           algorithm_params_t &algorithmParams)
{
    std::array<bool, CHAR_MAX> isAllowedElem{};
    for (const auto &alphabetElem : alphabet)
    {
        char elem = Rcpp::as<char>(alphabetElem);
        isAllowedElem[elem] = true;
    }

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult,
                         std::size_t seqBegin, std::size_t seqEnd) {
        KMerTaskConfig<EncodedStringVector> kMerTaskConfig(
            EncodedStringVector(isAllowedElem, sequences, seqBegin, seqEnd),
            config::DEFAULT_KMER_ITEM_SEPARATOR,
            config::DEFAULT_KMER_SECTION_SEPARATOR,
            userParams);
        updateKMerCountingResult<EncodedStringVector, kmer_dictionary_t>(
            kMerTaskConfig,
            algorithmParams,
            kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sequences.size(), userParams);
}

template <class algorithm_params_t,
          template <typename key, typename value, class...> class kmer_dictionary_t>
inline Rcpp::List parallelCountKMersSpecific(Rcpp::StringVector &sequences,
                                             Rcpp::StringVector &alphabet,
                                             const UserParams &userParams,
                                             algorithm_params_t &algorithmParams)
{
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
}

template <class algorithm_params_t,
          template <typename key, typename value, class...> class kmer_dictionary_t>
inline Rcpp::List sequentialCountKMersSpecific(Rcpp::StringVector &sequences,
                                               Rcpp::StringVector &alphabet,
                                               const UserParams &userParams,
                                               algorithm_params_t &algorithmParams)
{
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
}
