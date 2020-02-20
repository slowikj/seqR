#ifndef FIRST_UNORDERED_MAP_DICTIONARY_H
#define FIRST_UNORDERED_MAP_DICTIONARY_H

#include "dictionary.h"

#include <unordered_map>

template<class K, class V>
class unordered_map_dictionary : public dictionary<K, V> {
public:
  V &operator[](const K &key) override;
  
  [[nodiscard]] std::vector<int> get_keys() const override;
  
  unordered_map_dictionary() = default;
  
  unordered_map_dictionary(const unordered_map_dictionary<K, V> &) = default;
  
  unordered_map_dictionary<K, V> &operator=(const unordered_map_dictionary<K, V> &) = default;
  
  unordered_map_dictionary(unordered_map_dictionary<K, V> &&) noexcept = default;
  
  unordered_map_dictionary<K, V> &operator=(unordered_map_dictionary<K, V> &&) noexcept = default;
  
  ~unordered_map_dictionary() = default;
  
  bool is_present(const K &key) const override;
  
private:
  std::unordered_map<K, V> inner_map;
  
};

#include<algorithm>

template<class K, class V>
V &unordered_map_dictionary<K, V>::operator[](const K &key) {
  return this->inner_map[key];
}

template<class K, class V>
std::vector<int> unordered_map_dictionary<K, V>::get_keys() const {
  std::vector<int> res;
  res.reserve(this->inner_map.size());
  for (const std::pair<K, V> &elem: this->inner_map) {
    res.push_back(elem.first);
  }
  return res;
}

template<class K, class V>
bool unordered_map_dictionary<K, V>::is_present(const K &key) const {
  return this->inner_map.find(key) != std::end(this->inner_map);
}

#endif //FIRST_UNORDERED_MAP_DICTIONARY_H
