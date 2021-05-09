#pragma once

#include <memory>

#include "../hash/complex_hasher.h"
#include "../utils.h"
#include "globals.h"

namespace hashing {
template <class encoded_sequence_t>
class RollingWindow {
 public:
  using encoded_elem_t = typename encoded_sequence_t::encoded_elem_t;
  using hash_t = ComplexHasher::hash_t;

  RollingWindow(const encoded_sequence_t &sequence,
                ComplexHasher &&hasher) : sequence(sequence),
                                          hasher(std::move(hasher)) {
    this->nextElementIndex = 0;
  }

  RollingWindow() = delete;

  inline void resetIndex(int nextElementIndex) {
    this->nextElementIndex = nextElementIndex;
    util::clear<encoded_elem_t>(this->window);
    this->hasher.clear();
  }

  inline void append() {
    encoded_elem_t encodedElem = this->sequence[this->nextElementIndex];
    this->window.push(encodedElem);
    this->hasher.append(encodedElem);
    ++this->nextElementIndex;
  }

  inline void moveWindowRight() {
    removeFirst();
    append();
  }

  inline void removeFirst() {
    this->hasher.removeFirst(this->window.front());
    this->window.pop();
  }

  inline std::size_t sequenceSize() const {
    return this->sequence.size();
  }

  inline hash_t getWindowedHashes() const {
    return this->hasher.getHashes();
  }

  inline hash_t getWindowedPositionedHashes() const {
    return this->hasher.getHashes(this->nextElementIndex);
  }

  inline int currentBeginIndex() const {
    return this->nextElementIndex - this->window.size();
  }

 private:
  const encoded_sequence_t &sequence;

  ComplexHasher hasher;

  std::queue<encoded_elem_t> window;

  int nextElementIndex;
};
}  // namespace hashing
