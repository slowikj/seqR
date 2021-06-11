#pragma once

#include <Rcpp.h>

#include <array>
#include <limits>
#include <string>
#include <vector>

#include "../common_config.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"
#include "./encoded_sequence/raw_encoded_sequences_list.h"

template <class encoded_elem_t>
inline void addNewElement(
    const std::string &element,
    std::vector<encoded_elem_t> &encodedItems,
    std::unordered_map<std::string, encoded_elem_t> &alphabetEncoder,
    std::vector<std::string> &alphabetDecoder,
    encoded_elem_t invalidElemCode,
    encoded_elem_t &lastCnt,
    bool allElementsAllowed) {
  if (alphabetEncoder.find(element) != alphabetEncoder.end()) {
    encodedItems.push_back(alphabetEncoder[element]);
  } else {
    if (allElementsAllowed) {
      alphabetEncoder[element] = ++lastCnt;
      alphabetDecoder.push_back(element);
      encodedItems.push_back(lastCnt);
    } else {
      encodedItems.push_back(invalidElemCode);
    }
  }
}

template <class encoded_elem_t>
inline RawEncodedSequencesList<std::string, encoded_elem_t> encode(
    Rcpp::List sequences,
    std::size_t seqBegin,
    std::size_t seqEnd,
    std::unordered_map<std::string, encoded_elem_t> &alphabetEncoder,
    std::vector<std::string> &alphabetDecoder,
    encoded_elem_t invalidElemCode,
    encoded_elem_t &lastCnt,
    bool allElementsAllowed) {
  std::vector<encoded_elem_t> encodedItems{};
  std::vector<std::size_t> seqStarts{0};

  for (std::size_t i = seqBegin; i < seqEnd; ++i) {
    Rcpp::StringVector seq = sequences[i];
    for (const auto &seqElem : seq) {
      std::string cppElem = Rcpp::as<std::string>(seqElem);
      addNewElement(cppElem,
                    encodedItems,
                    alphabetEncoder, alphabetDecoder,
                    invalidElemCode, lastCnt, allElementsAllowed);
    }
    seqStarts.push_back(seqStarts.back() + seq.size());
  }

  return RawEncodedSequencesList<std::string, encoded_elem_t>(
      std::move(encodedItems),
      std::move(seqStarts),
      alphabetDecoder,
      invalidElemCode,
      allElementsAllowed);
}

template <class algorithm_params_t,
          class kmer_manager_t,
          template <typename key, typename value, class...> class result_dictionary_t>
inline Rcpp::List countKMersSpecific(Rcpp::List &sequences,
                                     Rcpp::StringVector &kmerAlphabet,
                                     const UserParams &userParams,
                                     algorithm_params_t &algorithmParams) {
  using encodedElemType = uint8_t;
  std::unordered_map<std::string, encodedElemType> alphabetEncoder{};
  std::vector<std::string> alphabetDecoder{"", ""};
  encodedElemType invalidElemCode = 1;
  encodedElemType encodingCnt = 1;
  bool allElementsAllowed = (kmerAlphabet[0] == config::ALPHABET_ALL_LABEL);
  if (!allElementsAllowed) {
    for (const auto &elem : kmerAlphabet) {
      std::string cppElem = Rcpp::as<std::string>(elem);
      alphabetEncoder[cppElem] = ++encodingCnt;
      alphabetDecoder.push_back(cppElem);
    }
  }

  auto batchFunc = [&](KMerCountingResult<result_dictionary_t> &kMerCountingResult,
                       std::size_t seqBegin, std::size_t seqEnd) {
    KMerTaskConfig<RawEncodedSequencesList<std::string, encodedElemType>> kMerTaskConfig(
        encode<encodedElemType>(sequences, seqBegin, seqEnd,
                                alphabetEncoder, alphabetDecoder,
                                invalidElemCode, encodingCnt,
                                allElementsAllowed),
        config::DEFAULT_KMER_ITEM_SEPARATOR,
        config::DEFAULT_KMER_SECTION_SEPARATOR,
        userParams);
    updateKMerCountingResult<RawEncodedSequencesList<std::string, encodedElemType>,
                             kmer_manager_t,
                             result_dictionary_t>(
        kMerTaskConfig,
        algorithmParams,
        kMerCountingResult);
  };

  return computeKMersInBatches<result_dictionary_t>(batchFunc, sequences.size(), userParams);
}
