#ifndef FIRST_UNORDEREDMAPDICTIONARY_H
#define FIRST_UNORDEREDMAPDICTIONARY_H

#include "dictionary.h"

#include <unordered_map>

template<class K, class V>
class UnorderedMapDictionary : public Dictionary<K, V> {
public:
  V &operator[](const K &key) override;
  
  [[nodiscard]] std::vector<int> get_keys() const override;
  
  UnorderedMapDictionary() = default;
  
  UnorderedMapDictionary(const UnorderedMapDictionary<K, V> &) = default;
  
  UnorderedMapDictionary<K, V> &operator=(const UnorderedMapDictionary<K, V> &) = default;
  
  UnorderedMapDictionary(UnorderedMapDictionary<K, V> &&) noexcept = default;
  
  UnorderedMapDictionary<K, V> &operator=(UnorderedMapDictionary<K, V> &&) noexcept = default;
  
  ~UnorderedMapDictionary() = default;
  
  bool is_present(const K &key) const override;
  
private:
  std::unordered_map<K, V> inner_map;
  
};

#include<algorithm>

template<class K, class V>
V &UnorderedMapDictionary<K, V>::operator[](const K &key) {
  return this->inner_map[key];
}

template<class K, class V>
std::vector<int> UnorderedMapDictionary<K, V>::get_keys() const {
  std::vector<int> res;
  res.reserve(this->inner_map.size());
  for (const std::pair<K, V> &elem: this->inner_map) {
    res.push_back(elem.first);
  }
  return res;
}

template<class K, class V>
bool UnorderedMapDictionary<K, V>::is_present(const K &key) const {
  return this->inner_map.find(key) != std::end(this->inner_map);
}

#endif //FIRST_UNORDEREDMAPDICTIONARY_H
