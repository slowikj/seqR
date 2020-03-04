#ifndef CUSTOM_DICTIONARY_H
#define CUSTOM_DICTIONARY_H

#include <vector>
#include <iterator>
#include <memory>
#include <iterator>
#include <unordered_map>
#include <functional>

template<class K, class V, class Hash=std::hash<K>>
class Dictionary {
public:
  
  class iterator {
  public:
    // iterator_traits
    using value_type = std::pair<const K&, V>;
    using pointer = std::pair<const K&, V> *;
    using reference = std::pair<const K&, V> &;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = void;
    
    explicit iterator(typename std::unordered_map<K, V, Hash>::iterator it) : container_iterator_(it) {}
    
    iterator& operator++() {
      ++this->container_iterator_;
      return *this;
    }
    
    iterator operator++(int) {
      auto old = iterator(this->container_iterator_);
      ++this->container_iterator_;
      return old;
    }
    
    value_type operator*() {
      return *(this->container_iterator_);
    }
    
    iterator operator->() {
      return (this->container_iterator_);
    }
    
    bool operator==(const iterator& other) {
      return this->container_iterator_ == other.container_iterator_;
    }
    
    bool operator!=(const iterator& other) {
      return !((*this) == other);
    }
    
  private:
    typename std::unordered_map<K, V, Hash>::iterator container_iterator_;
    
  };
  
  V &operator[](const K &key);
  
  std::vector<K> getKeys() const;
  
  std::size_t size() const;
  
  bool isPresent(const K &key) const;
  
  iterator begin();
  
  iterator end();
  
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
  typename Dictionary<K, V, Hash>::iterator Dictionary<K, V, Hash>::end() {
    return iterator(this->inner_map_.end());
  }

template<class K, class V, class Hash>
inline
  std::size_t Dictionary<K, V, Hash>::size() const {
    return this->inner_map_.size();
  }

#endif //CUSTOM_DICTIONARY_H
