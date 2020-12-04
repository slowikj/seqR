#ifndef SEQR_KMER_TASK_CONFIG_H
#define SEQR_KMER_TASK_CONFIG_H

#include <utility>
#include <vector>
#include "user_params.h"

template<class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

template<class input_elem_t>
using InputToStringItemConverter_t = std::function<std::string(const input_elem_t &)>;

template<class vector_t, class elem_t>
struct KMerTaskConfig {
    int sequencesNum;
    SequenceGetter_t<vector_t> sequenceGetter;
    InputToStringItemConverter_t<elem_t> inputToStringItemConverter;

    std::string kMerItemSeparator;
    std::string kMerSectionSeparator;

    const UserParams &userParams;

    KMerTaskConfig(int sequencesNum,
                   SequenceGetter_t<vector_t> sequenceGetter,
                   InputToStringItemConverter_t<elem_t> inputToStringItemConverter,
                   std::string kmerItemSeparator,
                   std::string kmerSectionSeparator,
                   const UserParams &userParams) :
            sequencesNum(sequencesNum),
            sequenceGetter(sequenceGetter),
            inputToStringItemConverter(inputToStringItemConverter),
            kMerItemSeparator(std::move(kmerItemSeparator)),
            kMerSectionSeparator(std::move(kmerSectionSeparator)),
            userParams(userParams) {
    }

};

#endif //SEQR_KMER_TASK_CONFIG_H
