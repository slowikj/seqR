#pragma once

#include <utility>
#include <vector>
#include "user_params.h"

template <class encoded_sequences_list_t>
struct KMerTaskConfig
{
    encoded_sequences_list_t EncodedSequencesList;

    std::string kMerItemSeparator;
    std::string kMerSectionSeparator;

    const UserParams &userParams;

    KMerTaskConfig(
        encoded_sequences_list_t &&EncodedSequencesList,
        std::string kmerItemSeparator,
        std::string kmerSectionSeparator,
        const UserParams &userParams)
        : EncodedSequencesList(std::move(EncodedSequencesList)),
          kMerItemSeparator(kmerItemSeparator),
          kMerSectionSeparator(kmerSectionSeparator),
          userParams(userParams)
    {
    }
};
