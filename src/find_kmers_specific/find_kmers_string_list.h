//
// Created by slowik on 28.11.2020.
//

#ifndef SEQR_FIND_KMERS_STRING_LIST_H
#define SEQR_FIND_KMERS_STRING_LIST_H

#include <Rcpp.h>
#include <vector>
#include "../kmer_task_config.h"
#include "../alphabet_encoder.h"
#include "../dictionary/unordered_map_wrapper.h"
#include "../kmer_counting_result.h"
#include "../kmer_task_solver.h"


template<class algorithm_params_t>
inline
Rcpp::List findKMersSpecific(Rcpp::List &sequences,
                             Rcpp::StringVector &alphabet,
                             std::vector<int> &gaps,
                             bool positionalKMers,
                             bool withKMerCounts,
                             const std::string &kmerDictionaryName,
                             int batchSize,
                             bool verbose,
                             algorithm_params_t &algorithmParams) {
    std::string alphabetStr("", alphabet.size());
    for (int i = 0; i < alphabet.size(); ++i) {
        alphabetStr[i] = Rcpp::as<char>(alphabet[i]);
    }
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<std::string, char, short, UnorderedMapWrapper>(alphabetStr));

    auto batchFunc = [&](KMerCountingResult &kMerCountingResult, int seqBegin, int seqEnd) {
        SafeSequencesStringListWrapper sequenceWrapper(sequences, seqBegin, seqEnd);
        KMerTaskConfig<SafeSequencesStringListWrapper::Row, char> kMerTaskConfig(
                (seqEnd - seqBegin),
                getStringSequenceGetter(sequenceWrapper),
                gaps,
                positionalKMers,
                withKMerCounts,
                getCharToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                SafeSequencesStringListWrapper::Row,
                char,
                short,
                AlphabetEncoding<char, short, UnorderedMapWrapper>,
                algorithm_params_t>(kMerTaskConfig,
                                    alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    };

    return computeKMersInBatches(batchFunc, sequences.size(), batchSize, verbose);
}


#endif //SEQR_FIND_KMERS_STRING_LIST_H
