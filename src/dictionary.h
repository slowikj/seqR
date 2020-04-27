#ifndef CUSTOM_DICTIONARY_H
#define CUSTOM_DICTIONARY_H

#include <vector>
#include <iterator>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <functional>
#include <algorithm>

template<class ITER_V, class INNER_CONTAINER_IT_T>
class iterator_t {
public:
    // iterator_traits
    using value_type = ITER_V;
    using pointer = ITER_V *;
    using reference = ITER_V &;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = void;

    explicit iterator_t(const INNER_CONTAINER_IT_T &it) : container_iterator_(it) {}

    inline iterator_t<ITER_V, INNER_CONTAINER_IT_T> &operator++() {
        ++this->container_iterator_;
        return *this;
    }

    inline iterator_t<ITER_V, INNER_CONTAINER_IT_T> operator++(int) {
        auto old = iterator(this->container_iterator_);
        ++this->container_iterator_;
        return old;
    }

    inline value_type operator*() {
        return *(this->container_iterator_);
    }

    inline iterator_t<ITER_V, INNER_CONTAINER_IT_T> operator->() {
        return (this->container_iterator_);
    }

    inline bool operator==(const iterator_t<ITER_V, INNER_CONTAINER_IT_T> &other) const {
        return this->container_iterator_ == other.container_iterator_;
    }

    inline bool operator!=(const iterator_t<ITER_V, INNER_CONTAINER_IT_T> &other) const {
        return !((*this) == other);
    }

private:
    INNER_CONTAINER_IT_T container_iterator_;

};

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
