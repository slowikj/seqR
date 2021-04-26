#pragma once

#include <utility>
#include <vector>
#include "user_params.h"

template <class encoded_sequences_t>
struct KMerTaskConfig
{
    encoded_sequences_t encodedSequences;

    std::string kMerItemSeparator;
    std::string kMerSectionSeparator;

    const UserParams &userParams;

    KMerTaskConfig(
        encoded_sequences_t &&encodedSequences,
        std::string kmerItemSeparator,
        std::string kmerSectionSeparator,
        const UserParams &userParams)
        : encodedSequences(std::move(encodedSequences)),
          kMerItemSeparator(kmerItemSeparator),
          kMerSectionSeparator(kmerSectionSeparator),
          userParams(userParams)
    {
    }
};
