#ifndef FIRST_DICTIONARY_H
#define FIRST_DICTIONARY_H

#include<vector>

template<class K, class V>
class dictionary {
public:
  virtual V &operator[](const K &key) = 0;
  
  [[nodiscard]] virtual std::vector<int> get_keys() const = 0;
  
  virtual bool is_present(const K &key) const = 0;
};


#endif //FIRST_DICTIONARY_H
