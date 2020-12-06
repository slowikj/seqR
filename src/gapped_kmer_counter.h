#ifndef GAPPED_KMER_COUNTER_H
#define GAPPED_KMER_COUNTER_H

#include "hash/polynomial_single_hasher.h"
#include "hash/types.h"
#include "kmer_manager.h"
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
    class PrefixSequencePolynomialHasher;

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
    inline hashing::multidim_hash_t getGappedKMerHashNotPositional(
            int beginPosition,
            const PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> &seqHasher,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals);

    template<class input_vector_t, class alphabet_encoding_t>
    inline hashing::multidim_hash_t getGappedKMerHash(
            int beginPosition,
            const PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> &seqHasher,
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
        PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> sequenceHasher(
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
                auto hash = std::move(getGappedKMerHash(seqInd, sequenceHasher, contiguousIntervals, isPositionalKMer));
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

    template<class input_vector_t,
            class alphabet_encoding_t>
    class PrefixSequencePolynomialHasher {
    public:
        PrefixSequencePolynomialHasher(input_vector_t &sequence,
                                       alphabet_encoding_t &alphabetEncoding,
                                       const std::vector<hashing::PolynomialSingleHasherConfig> &polynomialHasherConfigs)
                :
                polynomialHasherConfigs(polynomialHasherConfigs) {
            computePrefixValues(sequence, alphabetEncoding);
        }

        inline hashing::multidim_hash_t getHash(int begin, int end) const {
            hashing::multidim_hash_t res(polynomialHasherConfigs.size());
            for (int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
                int M = polynomialHasherConfigs[hasherInd].M;
                res[hasherInd] = ((prefixComplexHashes[end + 1][hasherInd]
                                   - (prefixComplexHashes[begin][hasherInd] *
                                      prefixP[end - begin + 1][hasherInd])) % M + M) % M;
            }
            return std::move(res);
        }

        inline int getHasherP(int hasherIndex, int power = 1) const {
            return this->prefixP[power][hasherIndex];
        }

        inline int getHasherM(int hasherIndex) const {
            return this->polynomialHasherConfigs[hasherIndex].M;
        }

        inline int getHashersNum() const {
            return this->polynomialHasherConfigs.size();
        }

    private:
        const std::vector<hashing::PolynomialSingleHasherConfig> &polynomialHasherConfigs;
        std::vector<std::vector<int>> prefixP;
        std::vector<hashing::multidim_hash_t> prefixComplexHashes;

        inline void computePrefixValues(input_vector_t &sequence,
                                        alphabet_encoding_t &alphabetEncoding) {
            initPrefixP(sequence.size(), polynomialHasherConfigs.size());
            initPrefixComplexHashes(sequence.size(), polynomialHasherConfigs.size());
            for (int seqInd = 0; seqInd < sequence.size(); ++seqInd) {
                appendPrefixValues(sequence, alphabetEncoding, seqInd);
            }
        }

        inline void initPrefixP(int sequenceLength, int hashersNum) {
            prefixP.reserve(sequenceLength);
            prefixP.push_back(std::move(std::vector<int>(hashersNum, 1)));
        }

        inline void initPrefixComplexHashes(int sequenceLength, int hashersNum) {
            prefixComplexHashes.reserve(sequenceLength);
            prefixComplexHashes.push_back(
                    std::move(hashing::multidim_hash_t(hashersNum))
            );
        }

        inline void appendPrefixValues(input_vector_t &sequence,
                                       alphabet_encoding_t &alphabetEncoding,
                                       int seqInd) {
            typename alphabet_encoding_t::encoded_elem_t encodedElem = alphabetEncoding.encode(sequence[seqInd]);
            appendCurrentComplexHash(encodedElem);
            appendCurrentPowerP();
        }

        inline void appendCurrentComplexHash(const typename alphabet_encoding_t::encoded_elem_t &encodedElem) {
            hashing::multidim_hash_t prefixHash(polynomialHasherConfigs.size());
            for (int hasherInd = 0; hasherInd < prefixHash.size(); ++hasherInd) {
                int P = polynomialHasherConfigs[hasherInd].P;
                int M = polynomialHasherConfigs[hasherInd].M;
                prefixHash[hasherInd] = ((prefixComplexHashes.back()[hasherInd]) * P + encodedElem) % M;
            }
            prefixComplexHashes.push_back(std::move(prefixHash));
        }

        inline void appendCurrentPowerP() {
            std::vector<int> powersP(polynomialHasherConfigs.size());
            for (int hasherInd = 0; hasherInd < powersP.size(); ++hasherInd) {
                int P = polynomialHasherConfigs[hasherInd].P;
                int M = polynomialHasherConfigs[hasherInd].M;
                powersP[hasherInd] = static_cast<int>(
                        (static_cast<long long>(prefixP.back()[hasherInd]) * P) % M);
            }
            prefixP.push_back(std::move(powersP));
        }
    };

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
        return std::move(res);
    }

    template<class input_vector_t, class alphabet_encoding_t>
    inline hashing::multidim_hash_t getGappedKMerHash(
            int beginPosition,
            const PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> &seqHasher,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
            bool isPositionalKMer) {
        hashing::multidim_hash_t res = std::move(
                getGappedKMerHashNotPositional<input_vector_t, alphabet_encoding_t>(
                        beginPosition, seqHasher, contiguousKMerIntervals)
        );
        if (isPositionalKMer) {
            res.push_back(beginPosition);
        }
        return std::move(res);
    }

    template<class input_vector_t, class alphabet_encoding_t>
    inline
    hashing::multidim_hash_t getGappedKMerHashNotPositional(
            int beginPosition,
            const PrefixSequencePolynomialHasher<input_vector_t, alphabet_encoding_t> &seqHasher,
            const std::vector<std::pair<int, int>> &contiguousKMerIntervals) {
        hashing::multidim_hash_t res(seqHasher.getHashersNum());
        for (const auto &kmerInterval: contiguousKMerIntervals) {
            auto intervalHash = std::move(seqHasher.getHash(
                    kmerInterval.first + beginPosition,
                    kmerInterval.second + beginPosition));
            int intervalLength = util::getIntervalLength(kmerInterval);
            for (int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
                int powerP = seqHasher.getHasherP(hasherInd, intervalLength);
                int M = seqHasher.getHasherM(hasherInd);
                res[hasherInd] = (((res[hasherInd]) * powerP + intervalHash[hasherInd]) % M + M) % M;
            }
        }
        return std::move(res);
    }

}

#endif
