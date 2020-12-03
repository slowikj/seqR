#ifndef SEQR_KMER_TASK_CONFIG_H
#define SEQR_KMER_TASK_CONFIG_H

#include <utility>
#include <vector>

template<class input_vector_t>
using SequenceGetter_t = std::function<input_vector_t(int)>;

template<class input_elem_t>
using InputToStringItemConverter_t = std::function<std::string(const input_elem_t &)>;

template<class vector_t, class elem_t>
struct KMerTaskConfig {
    int sequencesNum;
    SequenceGetter_t<vector_t> sequenceGetter;
    int k;
    bool positionalKMers;
    bool withKMerCounts;
    bool parallelMode;
    InputToStringItemConverter_t<elem_t> inputToStringItemConverter;
    std::vector<int> gaps;
    std::string kMerItemSeparator;
    std::string kMerSectionSeparator;

    KMerTaskConfig(int sequencesNum,
                   SequenceGetter_t<vector_t> sequenceGetter,
                   std::vector<int> &gaps,
                   bool positionalKMers,
                   bool withKMerCounts,
                   bool parallelMode,
                   InputToStringItemConverter_t<elem_t> inputToStringItemConverter,
                   std::string kmerItemSeparator,
                   std::string kmerSectionSeparator) :
            sequencesNum(sequencesNum),
            sequenceGetter(sequenceGetter),
            k(static_cast<int>(gaps.size()) + 1),
            gaps(gaps),
            positionalKMers(positionalKMers),
            withKMerCounts(withKMerCounts),
            parallelMode(parallelMode),
            inputToStringItemConverter(inputToStringItemConverter),
            kMerItemSeparator(std::move(kmerItemSeparator)),
            kMerSectionSeparator(std::move(kmerSectionSeparator)) {
    }

};

#endif //SEQR_KMER_TASK_CONFIG_H
