#ifndef KMER_COUNTS_MANAGER_H
#define KMER_COUNTS_MANAGER_H

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
        template<typename key, typename value> class kmer_dictionary_t>
class KMerManager {
public:
    using dict_t = kmer_dictionary_t<std::vector<int>, KMerHashInfo>;

    KMerManager() = default;

    KMerManager(const KMerManager<kmer_dictionary_t> &) = default;

    KMerManager(KMerManager<kmer_dictionary_t> &&) noexcept = default;

    KMerManager<kmer_dictionary_t> &
    operator=(const KMerManager<kmer_dictionary_t> &) = default;

    KMerManager<kmer_dictionary_t> &
    operator=(KMerManager<kmer_dictionary_t> &&) noexcept = default;

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
