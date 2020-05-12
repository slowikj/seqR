#ifndef GAPPED_KMER_COUNTER_H
#define GAPPED_KMER_COUNTER_H

#include "alphabet_encoder.h"
#include "hash/polynomial_single_hasher.h"
#include "kmer_counting_common_algorithm.h"
#include "kmer_manager.h"
#include "sequence_getter.h"
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

template<class vector_t>
inline std::vector<std::pair<int, int>> getContiguousIntervals(const vector_t &gaps) {
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

inline std::size_t getTotalKMerSize(const std::vector<int> &gaps) {
    return getSum(gaps) + gaps.size() + 1; // Rcpp::sum(gaps + 1) + 1
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline
std::vector<int> prepareNotAllowedItemsPrefixCount(
        input_vector_t &sequence,
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding) {
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

template<class input_vector_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
class PrefixSequencePolynomialHasher {
public:
    PrefixSequencePolynomialHasher(input_vector_t &sequence,
                                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                   const std::vector<PolynomialSingleHasherConfig> &polynomialHasherConfigs) :
            polynomialHasherConfigs(polynomialHasherConfigs) {
        computePrefixValues(sequence, alphabetEncoding);
    }

    inline std::vector<int> getHash(int begin, int end) const {
        std::vector<int> res(polynomialHasherConfigs.size());
        for (int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
            int M = polynomialHasherConfigs[hasherInd].M;
            res[hasherInd] = static_cast<int>(((prefixComplexHashes[end + 1][hasherInd]
                                                - (static_cast<long long>(prefixComplexHashes[begin][hasherInd]) *
                                                   prefixP[end - begin + 1][hasherInd])) % M + M) % M);
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
    const std::vector<PolynomialSingleHasherConfig> &polynomialHasherConfigs;
    std::vector<std::vector<int>> prefixP;
    std::vector<std::vector<int>> prefixComplexHashes;

    inline void computePrefixValues(input_vector_t &sequence,
                                    AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding) {
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
                std::move(std::vector<int>(hashersNum))
        );
    }

    inline void appendPrefixValues(input_vector_t &sequence,
                                   AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                   int seqInd) {
        encoded_elem_t encodedElem = alphabetEncoding.isAllowed(sequence[seqInd]) ?
                                     alphabetEncoding.encode(sequence[seqInd]) :
                                     alphabetEncoding.getNotAllowedEncodingNum();

        appendCurrentComplexHash(encodedElem);
        appendCurrentPowerP();
    }

    inline void appendCurrentComplexHash(const encoded_elem_t &encodedElem) {
        std::vector<int> prefixHash(polynomialHasherConfigs.size());
        for (int hasherInd = 0; hasherInd < prefixHash.size(); ++hasherInd) {
            int P = polynomialHasherConfigs[hasherInd].P;
            int M = polynomialHasherConfigs[hasherInd].M;
            prefixHash[hasherInd] = static_cast<int>(
                    (static_cast<long long>(prefixComplexHashes.back()[hasherInd]) * P + encodedElem) % M);
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

inline bool isGappedKMerAllowed(int seqBegin,
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

inline int getIntervalLength(const std::pair<int, int> &interval) {
    return interval.second - interval.first + 1;
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline
std::vector<int> getGappedKMerHashNotPositional(
        int beginPosition,
        const PrefixSequencePolynomialHasher<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t> &seqHasher,
        const std::vector<std::pair<int, int>> &contiguousKMerIntervals) {
    std::vector<int> res(seqHasher.getHashersNum());
    for (const auto &kmerInterval: contiguousKMerIntervals) {
        auto intervalHash = std::move(seqHasher.getHash(
                kmerInterval.first + beginPosition,
                kmerInterval.second + beginPosition));
        int intervalLength = getIntervalLength(kmerInterval);
        for (int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
            int powerP = seqHasher.getHasherP(hasherInd, intervalLength);
            int M = seqHasher.getHasherM(hasherInd);
            res[hasherInd] = static_cast<int>(
                    (static_cast<long long>(res[hasherInd]) * powerP + intervalHash[hasherInd]) % M);
        }
    }
    return std::move(res);
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline
std::vector<int> getGappedKMerHash(
        int beginPosition,
        const PrefixSequencePolynomialHasher<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t> &seqHasher,
        const std::vector<std::pair<int, int>> &contiguousKMerIntervals,
        bool isPositionalKMer) {
    std::vector<int> res = std::move(
            getGappedKMerHashNotPositional<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    beginPosition, seqHasher, contiguousKMerIntervals)
    );
    if (isPositionalKMer) {
        for (int hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
            res[hasherInd] = static_cast<int>(
                    (static_cast<long long>(res[hasherInd]) * seqHasher.getHasherP(hasherInd)
                     + beginPosition) % seqHasher.getHasherM(hasherInd));
        }
    }
    return std::move(res);
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
KMerManager<kmer_dictionary_t> countGappedKMers(const std::vector<int> &gaps,
                                                std::size_t totalKMerSize,
                                                input_vector_t &sequence,
                                                AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                                bool isPositionalKMer,
                                                bool withKMerCounts,
                                                const std::vector<PolynomialSingleHasherConfig> &hasherConfigs) {
    std::vector<std::pair<int, int>> contiguousIntervals = getContiguousIntervals(gaps);
    PrefixSequencePolynomialHasher<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t> sequenceHasher(
            sequence, alphabetEncoding, hasherConfigs
    );

    std::vector<int> notAllowedItemsPrefixCount = std::move(
            prepareNotAllowedItemsPrefixCount<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
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

template<class input_vector_t, class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
std::vector<KMerManager<kmer_dictionary_t>> parallelComputeGappedKMersCounts(
        KMerTaskConfig<input_vector_t, input_elem_t>& kMerTaskConfig,
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
        const std::vector<PolynomialSingleHasherConfig> &hasherConfigs) {
    std::size_t totalKMerSize = getTotalKMerSize(kMerTaskConfig.gaps);
    return std::move(
            parallelComputeKMers<input_vector_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    kMerTaskConfig.sequencesNum,
                    [&kMerTaskConfig, &alphabetEncoding, &totalKMerSize, &hasherConfigs]
                            (input_vector_t &v) -> KMerManager<kmer_dictionary_t> {
                        return countGappedKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                                kMerTaskConfig.gaps,
                                totalKMerSize,
                                v,
                                alphabetEncoding,
                                kMerTaskConfig.positionalKMers,
                                kMerTaskConfig.withKMerCounts,
                                hasherConfigs
                        );
                    },
                    kMerTaskConfig.sequenceGetter
            ));
}

#endif
