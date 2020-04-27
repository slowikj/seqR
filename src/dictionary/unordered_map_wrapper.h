#ifndef CUSTOM_DICTIONARY_H
#define CUSTOM_DICTIONARY_H

#include <vector>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include "iterator.h"

template<class K, class V, class Hash=std::hash<K>>
class UnorderedMapWrapper {
public:

    using iterator = iterator_t<std::pair<const K &, V>, typename std::unordered_map<K, V, Hash>::iterator>;
    using const_iterator = iterator_t<const std::pair<const K &, V>, typename std::unordered_map<K, V, Hash>::const_iterator>;

    inline V &operator[](const K &key) {
        return this->inner_map_[key];
    }

    inline V &operator[](K &&key) {
        return this->inner_map_[std::move(key)];
    }

    inline std::size_t size() const {
        return this->inner_map_.size();
    }

    inline bool isPresent(const K &key) const {
        return this->inner_map_.find(key) != std::end(this->inner_map_);
    }

    inline iterator begin() {
        return iterator(this->inner_map_.begin());
    }

    inline const_iterator begin() const noexcept {
        return const_iterator(this->inner_map_.begin());
    }

    inline iterator end() {
        return iterator(this->inner_map_.end());
    }

    inline const_iterator end() const noexcept {
        return const_iterator(this->inner_map_.end());
    }

    UnorderedMapWrapper() = default;

    UnorderedMapWrapper(const UnorderedMapWrapper<K, V, Hash> &) = default;

    UnorderedMapWrapper<K, V, Hash> &operator=(const UnorderedMapWrapper<K, V, Hash> &) = default;

    UnorderedMapWrapper(UnorderedMapWrapper<K, V, Hash> &&) noexcept = default;

    UnorderedMapWrapper<K, V, Hash> &operator=(UnorderedMapWrapper<K, V, Hash> &&) noexcept = default;

    ~UnorderedMapWrapper() = default;

private:
    std::unordered_map<K, V, Hash> inner_map_;
};

#endif
