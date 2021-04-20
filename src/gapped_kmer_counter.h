#pragma once

#include "hash/polynomial_single_hasher.h"
#include "hash/globals.h"
#include "kmer_manager.h"
#include "hash/prefix_sequence_polynomial_hasher.h"
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

namespace gappedKMers {

    template<class input_vector_t,
            class alphabet_encoding_t,
            template<typename key, typename value, typename...> class kmer_dictionary_t>
    inline KMerManager<kmer_dictionary_t> count(
            const std::vector<int> &gaps,
            std::size_t totalKMerSize,
            input_vector_t &sequence,
            alphabet_encoding_t &alphabetEncoding,
            bool isPositionalKMer,
            bool withKMerCounts,
            const std::vector<hashing::PolynomialSingleHasherConfig> &hasherConfigs);

    template<class input_vector_t,
            class alphabet_encoding_t>
    inline std::vector<int> prepareNotAllowedItemsPrefixCount(
            input_vector_t &sequence,
            alphabet_encoding_t &alphabetEncoding);

    template<class vector_t>
    inline std::vector<std::pair<int, int>> getContiguousIntervals(
            const vector_t &gaps);

    inline bool isGappedKMerAllowed(
            int seqBegin,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
            const std::vector<int> &notAllowedItemsPrefixCount);

    template<class input_vector_t, class alphabet_encoding_t>
    inline typename hashing::PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t>::hash_t
    getGappedKMerHash(
            int beginPosition,
            const hashing::PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> &seqHasher,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
            bool isPositionalKMer);

    // ------------------ IMPLEMENTATION ------------------

    template<class input_vector_t,
            class alphabet_encoding_t,
            template<typename key, typename value, typename...> class kmer_dictionary_t>
    inline
    KMerManager<kmer_dictionary_t> count(const std::vector<int> &gaps,
                                         std::size_t totalKMerSize,
                                         input_vector_t &sequence,
                                         alphabet_encoding_t &alphabetEncoding,
                                         bool isPositionalKMer,
                                         bool withKMerCounts,
                                         const std::vector<hashing::PolynomialSingleHasherConfig> &hasherConfigs) {
        std::vector<std::pair<int, int>> contiguousIntervals = getContiguousIntervals(gaps);
        hashing::PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> sequenceHasher(
                sequence, alphabetEncoding, hasherConfigs
        );

        std::vector<int> notAllowedItemsPrefixCount = std::move(
                prepareNotAllowedItemsPrefixCount<input_vector_t, alphabet_encoding_t>(
                        sequence,
                        alphabetEncoding
                )
        );

        KMerManager<kmer_dictionary_t> kMerManager(withKMerCounts);
        int limitSequenceIndex = static_cast<int>(sequence.size()) - totalKMerSize + 1;
        for (int seqInd = 0; seqInd < limitSequenceIndex; ++seqInd) {
            if (isGappedKMerAllowed(seqInd, contiguousIntervals, notAllowedItemsPrefixCount)) {
                auto hash = getGappedKMerHash(seqInd, sequenceHasher, contiguousIntervals, isPositionalKMer);
                kMerManager.add(std::move(hash), seqInd);
            }
        }

        return std::move(kMerManager);
    }

    template<class vector_t>
    inline std::vector<std::pair<int, int>> getContiguousIntervals(
            const vector_t &gaps) {
        std::vector<std::pair<int, int>> res;
        int currentKMerIndex = 0;
        for (int gapIndex = 0; gapIndex < gaps.size(); ++gapIndex) {
            int beginGapIndex = gapIndex;
            while (gapIndex < gaps.size() && gaps[gapIndex] == 0) {
                ++gapIndex;
            }
            res.emplace_back(currentKMerIndex,
                             currentKMerIndex + (gapIndex - beginGapIndex));

            if (gapIndex < gaps.size()) {
                currentKMerIndex += gaps[gapIndex] + (gapIndex - beginGapIndex) + 1;
            }

            if (gapIndex == gaps.size() - 1) {
                res.emplace_back(currentKMerIndex, currentKMerIndex);
            }
        }
        return res;
    }

    inline bool isGappedKMerAllowed(
            int seqBegin,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
            const std::vector<int> &notAllowedItemsPrefixCount) {
        return std::all_of(
                std::begin(contiguousKMerIntervals),
                std::end(contiguousKMerIntervals),
                [&notAllowedItemsPrefixCount, &seqBegin](const std::pair<int, int> &interval) -> bool {
                    return (notAllowedItemsPrefixCount[seqBegin + interval.second]
                            - (seqBegin + interval.first == 0 ?
                               0 :
                               notAllowedItemsPrefixCount[seqBegin + interval.first - 1])) == 0;
                }
        );
    }

    template<class input_vector_t, class alphabet_encoding_t>
    inline std::vector<int> prepareNotAllowedItemsPrefixCount(
            input_vector_t &sequence, alphabet_encoding_t &alphabetEncoding) {
        std::vector<int> res;
        res.reserve(sequence.size());
        for (int i = 0; i < sequence.size(); ++i) {
            bool isNotPresent = !alphabetEncoding.isAllowed(sequence[i]);
            res[i] = (i == 0) ?
                     isNotPresent :
                     res[i - 1] + isNotPresent;
        }
        return res;
    }

    template<class input_vector_t, class alphabet_encoding_t>
    inline typename hashing::PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t>::hash_t
    getGappedKMerHash(
            int beginPosition,
            const hashing::PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> &seqHasher,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
            bool isPositionalKMer) {
        auto res = seqHasher.getHashForSeveralIntervals(beginPosition, contiguousKMerIntervals);
        if (isPositionalKMer) {
            res.push_back(beginPosition);
        }
        return res;
    }
}
