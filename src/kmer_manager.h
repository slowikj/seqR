#ifndef KMER_COUNTS_MANAGER_H
#define KMER_COUNTS_MANAGER_H

#include "hash/custom_vector_hasher.h"
#include "hash/globals.h"
#include <vector>
#include <memory>

struct KMerHashInfo {
    int cnt;

    int seqStartPosition;

    explicit KMerHashInfo(int seqStartPosition) :
            cnt(0), seqStartPosition(seqStartPosition) {
    }

    KMerHashInfo(int seqStartPosition, int cnt) :
            cnt(cnt), seqStartPosition(seqStartPosition) {
    }

    KMerHashInfo() = default;

    KMerHashInfo &operator=(const KMerHashInfo &) = default;

    ~KMerHashInfo() = default;
};

template<
        template<typename key, typename value, typename...> class kmer_dictionary_t>
class KMerManager {
public:
    using dict_t = kmer_dictionary_t<hashing::config::multidim_hash_t, KMerHashInfo, hashing::config::multidim_hasher_t>;

    explicit KMerManager(bool kmerWithCounts) :
            kMerWithCounts(kmerWithCounts) {
    }

    KMerManager(const KMerManager<kmer_dictionary_t> &) = default;

    KMerManager() = default;

    KMerManager(KMerManager<kmer_dictionary_t> &&) noexcept = default;

    KMerManager<kmer_dictionary_t> &
    operator=(const KMerManager<kmer_dictionary_t> &) = default;

    KMerManager<kmer_dictionary_t> &
    operator=(KMerManager<kmer_dictionary_t> &&) noexcept = default;

    inline void add(hashing::config::multidim_hash_t &&hash, int position) {
        if (kMerWithCounts) {
            handleWithCounts(std::move(hash), position);
        } else {
            handleWithoutCounts(std::move(hash), position);
        }
    }

    inline const dict_t &getDictionary() const {
        return this->dictionary;
    }

private:
    dict_t dictionary;
    bool kMerWithCounts;

    inline void handleWithCounts(hashing::config::multidim_hash_t &&hash, int position) {
        if (!dictionary.isPresent(hash)) {
            dictionary[std::move(hash)] = KMerHashInfo(position, 1);
        } else {
            dictionary[hash].cnt++;
        }
    }

    inline void handleWithoutCounts(hashing::config::multidim_hash_t &&hash, int position) {
        dictionary[std::move(hash)] = KMerHashInfo(position, 1);
    }

};

#endif // KMER_COUNTS_MANAGER_H
