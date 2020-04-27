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
#include "kmer_counts_manager.h"
#include "kmer_counting_common.h"
#include "sequence_getter.h"
#include <vector>
#include <memory>
#include <functional>

template<class input_vector_t, class input_elem_t, class encoded_elem_t, template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline void updateKMerCounts(
        RollingWindow<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t> &rollingWindow,
        KMerCountsManager<kmer_dictionary_t> &kmerCountsManager,
        bool isPositionalKMer) {
    kmerCountsManager.add(
            isPositionalKMer ? rollingWindow.getWindowedPositionedHashes()
                             : rollingWindow.getWindowedHashes(),
            rollingWindow.currentBeginIndex()
    );
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t, template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline void countKMersForContiguousSeq(
        int k,
        int begin,
        int end,
        RollingWindow<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t> &rollingWindow,
        KMerCountsManager<kmer_dictionary_t> &kmerCountsManager,
        bool isPositionalKMer) {
    rollingWindow.resetIndex(begin);
    for (int i = 0; i < k; ++i) {
        rollingWindow.append();
    }
    for (int beginPosition = begin; beginPosition < end - k + 1; ++beginPosition) {
        updateKMerCounts(rollingWindow, kmerCountsManager, isPositionalKMer);
        rollingWindow.moveWindowRight();
    }
    updateKMerCounts(rollingWindow, kmerCountsManager, isPositionalKMer);
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t, template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline std::vector<int> computeNotAllowedPositions(
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
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

template<class input_vector_t, class input_elem_t, class encoded_elem_t, template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline KMerCountsManager<kmer_dictionary_t> countKMers(
        int k,
        input_vector_t &sequence,
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
        bool isPositionalKMer,
        ComplexHasher &&complexHasher) {
    KMerCountsManager<kmer_dictionary_t> kmerCountsManager;
    RollingWindow<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t> rollingWindow(
            sequence, std::move(complexHasher), alphabetEncoding
    );
    auto notAllowedSequencePositions = computeNotAllowedPositions(alphabetEncoding, sequence);
    for (int i = 0; i < notAllowedSequencePositions.size() - 1; ++i) {
        int allowedItemsBetween = notAllowedSequencePositions[i + 1] - notAllowedSequencePositions[i] - 1;
        if (allowedItemsBetween >= k) {
            int begin = notAllowedSequencePositions[i] + 1;
            int end = notAllowedSequencePositions[i + 1] - 1;
            countKMersForContiguousSeq<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    k, begin, end, rollingWindow, kmerCountsManager, isPositionalKMer
            );
        }
    }
    return std::move(kmerCountsManager);
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t, template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
std::vector<KMerCountsManager<kmer_dictionary_t>> parallelComputeKMerCounts(
        int k,
        bool positionalKMer,
        int sequencesNum,
        SequenceGetter_t<input_vector_t> sequenceGetter,
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
        std::function<ComplexHasher()> complexHasherFactory) {
    return std::move(
            parallelComputeKMerCounts<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    sequencesNum,
                    [k, positionalKMer, &alphabetEncoding, &complexHasherFactory]
                            (input_vector_t &v) -> KMerCountsManager<kmer_dictionary_t> {
                        return countKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                                k,
                                v,
                                alphabetEncoding,
                                positionalKMer,
                                std::move(complexHasherFactory())
                        );
                    },
                    sequenceGetter
            ));
}

#endif //KMER_COUNTER_H
