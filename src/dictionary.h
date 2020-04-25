#ifndef CUSTOM_DICTIONARY_H
#define CUSTOM_DICTIONARY_H

#include <vector>
#include <iterator>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <functional>

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

    iterator_t<ITER_V, INNER_CONTAINER_IT_T> &operator++() {
        ++this->container_iterator_;
        return *this;
    }

    iterator_t<ITER_V, INNER_CONTAINER_IT_T> operator++(int) {
        auto old = iterator(this->container_iterator_);
        ++this->container_iterator_;
        return old;
    }

    value_type operator*() {
        return *(this->container_iterator_);
    }

    iterator_t<ITER_V, INNER_CONTAINER_IT_T> operator->() {
        return (this->container_iterator_);
    }

    bool operator==(const iterator_t<ITER_V, INNER_CONTAINER_IT_T> &other) const {
        return this->container_iterator_ == other.container_iterator_;
    }

    bool operator!=(const iterator_t<ITER_V, INNER_CONTAINER_IT_T> &other) const {
        return !((*this) == other);
    }

private:
    INNER_CONTAINER_IT_T container_iterator_;

};

template<class K, class V, class Hash=std::hash<K>>
class Dictionary {
public:

    using iterator = iterator_t<std::pair<const K &, V>, typename std::unordered_map<K, V, Hash>::iterator>;
    using const_iterator = iterator_t<const std::pair<const K &, V>, typename std::unordered_map<K, V, Hash>::const_iterator>;

    V &operator[](const K &key);

    V &operator[](K &&key);

    std::vector<K> getKeys() const;

    std::size_t size() const;

    bool isPresent(const K &key) const;

    iterator begin();

    const_iterator begin() const noexcept;

    iterator end();

    const_iterator end() const noexcept;

    Dictionary() = default;

    Dictionary(const Dictionary<K, V, Hash> &) = default;

    Dictionary<K, V, Hash> &operator=(const Dictionary<K, V, Hash> &) = default;

    Dictionary(Dictionary<K, V, Hash> &&) noexcept = default;

    Dictionary<K, V, Hash> &operator=(Dictionary<K, V, Hash> &&) noexcept = default;

    ~Dictionary() = default;

private:
    std::unordered_map<K, V, Hash> inner_map_;
};

#include <algorithm>

template<class K, class V, class Hash>
inline
V &Dictionary<K, V, Hash>::operator[](const K &key) {
    return this->inner_map_[key];
}

template<class K, class V, class Hash>
inline
V &Dictionary<K, V, Hash>::operator[](K &&key) {
    return this->inner_map_[std::move(key)];
}

template<class K, class V, class Hash>
inline
std::vector<K> Dictionary<K, V, Hash>::getKeys() const {
    std::vector<K> res;
    res.reserve(this->inner_map_.size());
    for (const std::pair<K, V> &elem: this->inner_map_) {
        res.push_back(elem.first);
    }
    return res;
}

template<class K, class V, class Hash>
inline
bool Dictionary<K, V, Hash>::isPresent(const K &key) const {
    return this->inner_map_.find(key) != std::end(this->inner_map_);
}

template<class K, class V, class Hash>
inline
typename Dictionary<K, V, Hash>::iterator Dictionary<K, V, Hash>::begin() {
    return iterator(this->inner_map_.begin());
}

template<class K, class V, class Hash>
inline
typename Dictionary<K, V, Hash>::const_iterator Dictionary<K, V, Hash>::begin() const noexcept {
    return const_iterator(this->inner_map_.begin());
}

template<class K, class V, class Hash>
inline
typename Dictionary<K, V, Hash>::iterator Dictionary<K, V, Hash>::end() {
    return iterator(this->inner_map_.end());
}

template<class K, class V, class Hash>
inline
typename Dictionary<K, V, Hash>::const_iterator Dictionary<K, V, Hash>::end() const noexcept {
    return const_iterator(this->inner_map_.end());
}

template<class K, class V, class Hash>
inline
std::size_t Dictionary<K, V, Hash>::size() const {
    return this->inner_map_.size();
}

#endif
