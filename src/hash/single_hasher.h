#ifndef SECOND_SINGLE_HASHER_H
#define SECOND_SINGLE_HASHER_H

#include<queue>
#include<functional>
#include<utility>

class SingleHasher {
public:
  
  virtual void append(const int &elem) = 0;
  
  virtual void removeFirst(const int &item) = 0;
  
  [[nodiscard]] virtual int getHash() const {
    return currentHash;
  }
  
  [[nodiscard]] virtual int getHash(int position) const = 0;
  
protected:
  int currentHash = 0;
  
};

#endif //SECOND_SINGLE_HASHER_H
