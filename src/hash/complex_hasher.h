#ifndef SECOND_COMPLEX_HASHER_H
#define SECOND_COMPLEX_HASHER_H

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>
#include "single_hasher.h"

class ComplexHasher {
public:
  
  explicit ComplexHasher(std::vector<std::unique_ptr<SingleHasher>> &&singleHashers) :
    singleHashers(std::move(singleHashers)) {
  }
  
  void clear() {
    std::for_each(std::begin(this->singleHashers), std::end(this->singleHashers),
                  [](std::unique_ptr<SingleHasher> &singleHasher) {
                    singleHasher->clear();
                  }
    );
  }
  
  void append(const int &elem) {
    std::for_each(std::begin(this->singleHashers), std::end(this->singleHashers),
                  [&elem](std::unique_ptr<SingleHasher> &singleHasher) {
                    singleHasher->append(elem);
                  });
  }
  
  void removeFirst(const int &elem) {
    std::for_each(std::begin(this->singleHashers), std::end(this->singleHashers),
                  [&elem](std::unique_ptr<SingleHasher> &singleHasher) {
                    singleHasher->removeFirst(elem);
                  });
  }
  
  std::vector<int> getHashes(int position) const {
    return prepareResultHashes(
      [&position](const std::unique_ptr<SingleHasher> &singleHasher) -> int {
        return singleHasher->getHash(position);
      }
    );
  }
  
  std::vector<int> getHashes() const {
    return prepareResultHashes(
      [](const std::unique_ptr<SingleHasher> &singleHasher) -> int {
        return singleHasher->getHash();
      }
    );
  }
  
private:
  std::vector<std::unique_ptr<SingleHasher>> singleHashers;
  
  std::vector<int> prepareResultHashes(
      std::function<int(const std::unique_ptr<SingleHasher> &)> &&transformFunc) const {
    std::vector<int> resultHashes;
    std::transform(std::begin(singleHashers), std::end(singleHashers),
                   std::back_inserter(resultHashes),
                   transformFunc);
    return resultHashes;
  }
  
};

#endif //SECOND_COMPLEX_HASHER_H
