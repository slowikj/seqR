#pragma once

#include <utility>
#include <vector>
#include "user_params.h"
#include "encoded_sequence.h"

template<class init_elem_t, class encoded_elem_t>
struct KMerTaskConfig {
    EncodedSequences<init_elem_t, encoded_elem_t> encodedSequences;

    std::string kMerItemSeparator;
    std::string kMerSectionSeparator;

    const UserParams &userParams;

    KMerTaskConfig(EncodedSequences<init_elem_t, encoded_elem_t> &&encodedSequences,
                   std::string kmerItemSeparator,
                   std::string kmerSectionSeparator,
                   const UserParams &userParams) :
            encodedSequences(std::move(encodedSequences)),
            kMerItemSeparator(kmerItemSeparator),
            kMerSectionSeparator(kmerSectionSeparator),
            userParams(userParams) {
    }

};
