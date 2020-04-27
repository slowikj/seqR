#ifndef KMER_COUNTS_MANAGER_H
#define KMER_COUNTS_MANAGER_H

#include "dictionary.h"
#include "hash/custom_hashers.h"
#include <vector>
#include <memory>

struct KMerHashInfo {

    int cnt;

    int seqStartPosition;

    explicit KMerHashInfo(int seqStartPosition) :
            cnt(0), seqStartPosition(seqStartPosition) {
    }

    KMerHashInfo(int seqStartPosition, int cnt) :
            seqStartPosition(seqStartPosition), cnt(cnt) {
    }

    KMerHashInfo() = default;

    KMerHashInfo &operator=(const KMerHashInfo &) = default;

    ~KMerHashInfo() = default;
};

template<
        template<typename key, typename value, typename...> class kmer_dictionary_t>
class KMerCountsManager {
public:
    using dict_t = kmer_dictionary_t<std::vector<int>, KMerHashInfo>;

    KMerCountsManager() = default;

    KMerCountsManager(const KMerCountsManager<kmer_dictionary_t> &) = default;

    KMerCountsManager(KMerCountsManager<kmer_dictionary_t> &&) noexcept = default;

    KMerCountsManager<kmer_dictionary_t> &
    operator=(const KMerCountsManager<kmer_dictionary_t> &) = default;

    KMerCountsManager<kmer_dictionary_t> &
    operator=(KMerCountsManager<kmer_dictionary_t> &&) noexcept = default;

    inline void add(std::vector<int> &&hash, int position) {
        if (!this->dictionary.isPresent(hash)) {
            this->dictionary[std::move(hash)] = KMerHashInfo(position, 1);
        } else {
            this->dictionary[hash].cnt++;
        }
    }

    inline const dict_t &getDictionary() const {
        return this->dictionary;
    }

private:
    dict_t dictionary;

};

#endif // KMER_COUNTS_MANAGER_H
