#pragma once

#include "hash/custom_vector_hasher.h"
#include "hash/globals.h"
#include <vector>
#include <memory>

struct KMerHashInfo
{
    int cnt;

    int seqStartPosition;

    explicit KMerHashInfo(int seqStartPosition) : cnt(0), seqStartPosition(seqStartPosition)
    {
    }

    KMerHashInfo(int seqStartPosition, int cnt) : cnt(cnt), seqStartPosition(seqStartPosition)
    {
    }

    KMerHashInfo() = default;

    KMerHashInfo &operator=(const KMerHashInfo &) = default;

    ~KMerHashInfo() = default;
};

template <
    template <typename key, typename value, typename...> class kmer_dictionary_t_>
class CountingKMerManager
{
public:
    using hash_t = hashing::config::multidim_hash_t;
    using dict_t = kmer_dictionary_t_<hash_t, KMerHashInfo, hashing::config::multidim_hasher_t>;

    CountingKMerManager(const CountingKMerManager &) = default;

    CountingKMerManager() = default;

    CountingKMerManager(CountingKMerManager &&) = default;

    CountingKMerManager &
    operator=(const CountingKMerManager &) = default;

    CountingKMerManager &
    operator=(CountingKMerManager &&) = default;

    inline void add(hash_t &&hash, int position)
    {
        if (!dictionary.isPresent(hash))
        {
            dictionary[std::move(hash)] = KMerHashInfo(position, 1);
        }
        else
        {
            dictionary[hash].cnt++;
        }
    }

    inline const dict_t &getDictionary() const
    {
        return this->dictionary;
    }

private:
    dict_t dictionary;
};

template <
    template <typename key, typename value, typename...> class kmer_dictionary_t_>
class PresenceKMerManager
{
public:
    using hash_t = hashing::config::multidim_hash_t;
    using dict_t = kmer_dictionary_t_<hash_t, KMerHashInfo, hashing::config::multidim_hasher_t>;

    PresenceKMerManager(const PresenceKMerManager &) = default;

    PresenceKMerManager() = default;

    PresenceKMerManager(PresenceKMerManager &&) = default;

    PresenceKMerManager &
    operator=(const PresenceKMerManager &) = default;

    PresenceKMerManager &
    operator=(PresenceKMerManager &&) = default;

    inline void add(hash_t &&hash, int position)
    {
        dictionary[std::move(hash)] = KMerHashInfo(position, 1);
    }

    inline const dict_t &getDictionary() const
    {
        return this->dictionary;
    }

private:
    dict_t dictionary;
};
