#pragma once

#include <utility>
#include <vector>

#include "user_params.h"

template <class encoded_sequences_list_t>
struct KMerTaskConfig {
  encoded_sequences_list_t encodedSequencesList;

  std::string kMerItemSeparator;
  std::string kMerSectionSeparator;

  const UserParams &userParams;

  KMerTaskConfig(
      encoded_sequences_list_t &&encodedSequencesList,
      std::string kmerItemSeparator,
      std::string kmerSectionSeparator,
      const UserParams &userParams)
      : encodedSequencesList(std::move(encodedSequencesList)),
        kMerItemSeparator(kmerItemSeparator),
        kMerSectionSeparator(kmerSectionSeparator),
        userParams(userParams) {}
};
