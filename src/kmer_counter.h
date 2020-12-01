#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include "hash/rolling_window.h"
#include "hash/complex_hasher.h"
#include "hash/single_hasher.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_manager.h"
#include "kmer_counting_common_algorithm.h"
#include "sequence_getter.h"
#include <vector>
#include <memory>
#include <functional>

template<class input_vector_t,
        class alphabet_encoding_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline void updateKMers(
        RollingWindow<input_vector_t, alphabet_encoding_t> &rollingWindow,
        KMerManager<kmer_dictionary_t> &kMerManager,
        bool isPositionalKMer) {
    kMerManager.add(
            isPositionalKMer ? rollingWindow.getWindowedPositionedHashes()
                             : rollingWindow.getWindowedHashes(),
            rollingWindow.currentBeginIndex()
    );
}

template<class input_vector_t,
        class alphabet_encoding_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline void countKMersForContiguousSeq(
        int k,
        int begin,
        int end,
        RollingWindow<input_vector_t, alphabet_encoding_t> &rollingWindow,
        KMerManager<kmer_dictionary_t> &kMerManager,
        bool isPositionalKMer) {
    rollingWindow.resetIndex(begin);
    for (int i = 0; i < k; ++i) {
        rollingWindow.append();
    }
    for (int beginPosition = begin; beginPosition < end - k + 1; ++beginPosition) {
        updateKMers(rollingWindow, kMerManager, isPositionalKMer);
        rollingWindow.moveWindowRight();
    }
    updateKMers(rollingWindow, kMerManager, isPositionalKMer);
}

template<class input_vector_t,
        class alphabet_encoding_t>
inline std::vector<int> computeNotAllowedPositions(
        alphabet_encoding_t &alphabetEncoding,
        input_vector_t &sequence) {
    std::vector<int> res;
    res.push_back(-1); // left sentinel
    for (int seq_i = 0; seq_i < sequence.size(); ++seq_i) {
        if (!alphabetEncoding.isAllowed(sequence[seq_i])) {
            res.push_back(seq_i);
        }
    }
    res.push_back(sequence.size()); // right sentinel
    return res;
}

template<class input_vector_t,
        class alphabet_encoding_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline KMerManager<kmer_dictionary_t> countKMers(
        int k,
        input_vector_t &sequence,
        alphabet_encoding_t &alphabetEncoding,
        bool isPositionalKMer,
        bool withKMerCounts,
        ComplexHasher &&complexHasher) {
    KMerManager<kmer_dictionary_t> kMerManager(withKMerCounts);
    RollingWindow<input_vector_t, alphabet_encoding_t> rollingWindow(
            sequence, std::move(complexHasher), alphabetEncoding
    );
    auto notAllowedSequencePositions = computeNotAllowedPositions(alphabetEncoding, sequence);
    for (int i = 0; i < notAllowedSequencePositions.size() - 1; ++i) {
        int allowedItemsBetween = notAllowedSequencePositions[i + 1] - notAllowedSequencePositions[i] - 1;
        if (allowedItemsBetween >= k) {
            int begin = notAllowedSequencePositions[i] + 1;
            int end = notAllowedSequencePositions[i + 1] - 1;
            countKMersForContiguousSeq<input_vector_t, alphabet_encoding_t, kmer_dictionary_t>(
                    k, begin, end, rollingWindow, kMerManager, isPositionalKMer
            );
        }
    }
    return std::move(kMerManager);
}

#endif //KMER_COUNTER_H
