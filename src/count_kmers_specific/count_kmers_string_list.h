#pragma once

#include <Rcpp.h>
#include <vector>
#include <limits>
#include <string>
#include <array>
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "../encoded_sequence/raw_encoded_sequences_list.h"
#include "../common_config.h"

template <class encoded_elem_t>
inline RawEncodedSequencesList<std::string, encoded_elem_t> encode(
    Rcpp::List sequences,
    int seqBegin,
    int seqEnd,
    std::unordered_map<std::string, encoded_elem_t> &alphabetEncoder,
    const std::vector<std::string> &alphabetDecoder,
    encoded_elem_t invalidElemCode = 1)
{
    std::vector<encoded_elem_t> encodedItems{};
    std::vector<std::size_t> seqStarts{0};

    for (std::size_t i = seqBegin; i < seqEnd; ++i)
    {
        Rcpp::StringVector seq = sequences[i];
        for (const auto &seqElem : seq)
        {
            std::string cppElem = Rcpp::as<std::string>(seqElem);
            if (alphabetEncoder.find(cppElem) != alphabetEncoder.end())
            {
                encodedItems.push_back(invalidElemCode);
            }
            else
            {
                encodedItems.push_back(alphabetEncoder[cppElem]);
            }
        }
        seqStarts.push_back(seqStarts.back() + seq.size());
    }

    return RawEncodedSequencesList<std::string, encoded_elem_t>(
        std::move(encodedItems),
        std::move(seqStarts),
        alphabetDecoder,
        invalidElemCode
    );
}

template <class algorithm_params_t,
          template <typename key, typename value, class...> class kmer_dictionary_t>
inline Rcpp::List commonCountKMersSpecific(Rcpp::List &sequences,
                                           Rcpp::StringVector &alphabet,
                                           const UserParams &userParams,
                                           algorithm_params_t &algorithmParams)
{
    using encodedElemType = uint8_t;
    std::unordered_map<std::string, encodedElemType> alphabetEncoder{};
    std::vector<std::string> alphabetDecoder{"", ""};
    encodedElemType cnt = 1;
    for (const auto &elem : alphabet)
    {
        std::string cppElem = Rcpp::as<std::string>(elem);
        alphabetEncoder[cppElem] = ++cnt;
        alphabetDecoder.push_back(cppElem);
    }

    auto batchFunc = [&](KMerCountingResult<kmer_dictionary_t> &kMerCountingResult,
                         int seqBegin, int seqEnd) {
        KMerTaskConfig<RawEncodedSequencesList<std::string, encodedElemType>> kMerTaskConfig(
            encode<encodedElemType>(sequences, seqBegin, seqEnd, alphabetEncoder, alphabetDecoder),
            config::DEFAULT_KMER_ITEM_SEPARATOR,
            config::DEFAULT_KMER_SECTION_SEPARATOR,
            userParams);
        updateKMerCountingResult<RawEncodedSequencesList<std::string, encodedElemType>,
                                 kmer_dictionary_t>(
            kMerTaskConfig,
            algorithmParams,
            kMerCountingResult);
    };

    return computeKMersInBatches<kmer_dictionary_t>(batchFunc, sequences.size(), userParams);
}

template <class algorithm_params_t,
          template <typename key, typename value, class...> class kmer_dictionary_t>
inline Rcpp::List parallelCountKMersSpecific(Rcpp::List &sequences,
                                             Rcpp::StringVector &alphabet,
                                             const UserParams &userParams,
                                             algorithm_params_t &algorithmParams)
{
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
}

template <class algorithm_params_t,
          template <typename key, typename value, class...> class kmer_dictionary_t>
inline Rcpp::List sequentialCountKMersSpecific(Rcpp::List &sequences,
                                               Rcpp::StringVector &alphabet,
                                               const UserParams &userParams,
                                               algorithm_params_t &algorithmParams)
{
    return commonCountKMersSpecific<algorithm_params_t, kmer_dictionary_t>(
        sequences, alphabet, userParams, algorithmParams);
}