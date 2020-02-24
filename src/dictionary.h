#ifndef CUSTOM_DICTIONARY_H
#define CUSTOM_DICTIONARY_H

#include <vector>
#include <iterator>
#include <memory>
#include <iterator>
#include <unordered_map>

template<class K, class V>
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
    
    explicit iterator(typename std::unordered_map<K, V>::iterator it) : container_iterator_(it) {}
    
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
    typename std::unordered_map<K, V>::iterator container_iterator_;
    
  };
  
  V &operator[](const K &key);
  
  const V &operator[](const K &key) const;
  
  std::vector<K> getKeys() const;
  
  std::size_t size() const;
  
  bool isPresent(const K &key) const;
  
  iterator begin();
  
  iterator end();
  
  Dictionary() = default;
  
  Dictionary(const Dictionary<K, V> &) = default;
  
  Dictionary<K, V> &operator=(const Dictionary<K, V> &) = default;
  
  Dictionary(Dictionary<K, V> &&) noexcept = default;
  
  Dictionary<K, V> &operator=(Dictionary<K, V> &&) noexcept = default;
  
  ~Dictionary() = default;
  
private:
  std::unordered_map<K, V> inner_map_;
};

#include <algorithm>

template<class K, class V>
inline
  V &Dictionary<K, V>::operator[](const K &key) {
    return this->inner_map_[key];
  }

template<class K, class V>
inline
  const V &Dictionary<K, V>::operator[](const K &key) const {
    return this->inner_map_[key];
  }

template<class K, class V>
inline
  std::vector<K> Dictionary<K, V>::getKeys() const {
    std::vector<K> res;
    res.reserve(this->inner_map_.size());
    for (const std::pair<K, V> &elem: this->inner_map_) {
      res.push_back(elem.first);
    }
    return res;
  }

template<class K, class V>
inline
  bool Dictionary<K, V>::isPresent(const K &key) const {
    return this->inner_map_.find(key) != std::end(this->inner_map_);
  }

template<class K, class V>
inline
  typename Dictionary<K, V>::iterator Dictionary<K, V>::begin() {
    return iterator(this->inner_map_.begin());
  }

template<class K, class V>
inline
  typename Dictionary<K, V>::iterator Dictionary<K, V>::end() {
    return iterator(this->inner_map_.end());
  }

template<class K, class V>
inline
  std::size_t Dictionary<K, V>::size() const {
    return this->inner_map_.size();
  }

#endif //CUSTOM_DICTIONARY_H
