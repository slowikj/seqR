#ifndef SEQR_KMER_TASK_CONFIG_H
#define SEQR_KMER_TASK_CONFIG_H

#include "sequence_getter.h"
#include "input_to_string_item_converter.h"
#include <utility>
#include <vector>

const std::string DEFAULT_KMER_ITEM_SEPARATOR = ".";

const std::string DEFAULT_KMER_SECTION_SEPARATOR = "_";

template<class vector_t, class elem_t>
struct KMerTaskConfig {
    int sequencesNum;
    SequenceGetter_t<vector_t> sequenceGetter;
    int k;
    bool positionalKMers;
    bool withKMerCounts;
    InputToStringItemConverter_t<elem_t> inputToStringItemConverter;
    std::vector<int> gaps;
    std::string kMerItemSeparator;
    std::string kMerSectionSeparator;

    KMerTaskConfig(int sequencesNum,
                   SequenceGetter_t<vector_t> sequenceGetter,
                   std::vector<int>& gaps,
                   bool positionalKMers,
                   bool withKMerCounts,
                   InputToStringItemConverter_t<elem_t> inputToStringItemConverter,
                   std::string kmerItemSeparator = DEFAULT_KMER_ITEM_SEPARATOR,
                   std::string kmerSectionSeparator = DEFAULT_KMER_SECTION_SEPARATOR):
            sequencesNum(sequencesNum),
            sequenceGetter(sequenceGetter),
            k(static_cast<int>(gaps.size()) + 1),
            gaps(gaps),
            positionalKMers(positionalKMers),
            withKMerCounts(withKMerCounts),
            inputToStringItemConverter(inputToStringItemConverter),
            kMerItemSeparator(std::move(kmerItemSeparator)),
            kMerSectionSeparator(std::move(kmerSectionSeparator)) {
    }

    KMerTaskConfig(int sequencesNum,
                   SequenceGetter_t<vector_t> sequenceGetter,
                   int k,
                   bool positionalKMers,
                   bool withKMerCounts,
                   InputToStringItemConverter_t<elem_t> inputToStringItemConverter,
                   std::string kmerItemSeparator = DEFAULT_KMER_ITEM_SEPARATOR,
                   std::string kmerSectionSeparator = DEFAULT_KMER_SECTION_SEPARATOR) :
            sequencesNum(sequencesNum),
            sequenceGetter(sequenceGetter),
            k(k),
            positionalKMers(positionalKMers),
            withKMerCounts(withKMerCounts),
            inputToStringItemConverter(inputToStringItemConverter),
            gaps(k - 1),
            kMerItemSeparator(std::move(kmerItemSeparator)),
            kMerSectionSeparator(std::move(kmerSectionSeparator)) {
    }
};

#endif //SEQR_KMER_TASK_CONFIG_H
