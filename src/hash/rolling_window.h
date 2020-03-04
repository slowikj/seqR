#ifndef ROLLING_WINDOW_H
#define ROLLING_WINDOW_H

#include <memory>
#include "../alphabet_encoder/alphabet_encoder.h"
#include "../hash/complex_hasher.h"
#include "../utils.h"

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
class RollingWindow {
public:
  
  RollingWindow(input_vector_t& sequence,
                ComplexHasher&& hasher,
                AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding):
    sequence(sequence),
    hasher(std::move(hasher)),
    alphabetEncoding(alphabetEncoding) {
    this->nextElementIndex = 0;
  }
  
  RollingWindow() = delete;
  
  void resetIndex(int nextElementIndex) {
    this->nextElementIndex = nextElementIndex;
    clear(this->window);
    this->hasher.clear();
  }
  
  void append() {
    encoded_elem_t encodedElem = this->alphabetEncoding.encode(
      this->sequence[this->nextElementIndex]
    );
    this->window.push(encodedElem);
    this->hasher.append(encodedElem);
    ++this->nextElementIndex;
  }
  
  void moveWindowRight() {
    removeFirst();
    append();
  }
  
  void removeFirst() {
    this->hasher.removeFirst(this->window.front());
    this->window.pop();
  }
  
  std::size_t sequenceSize() const {
    return this->sequence.size();
  }
  
  std::vector<int> getWindowedHashes() const {
    return this->hasher.getHashes();
  }
  
  std::vector<int> getWindowedPositionedHashes() const {
    return this->hasher.getHashes(this->nextElementIndex);
  }
  
  int currentBeginIndex() const {
    return this->nextElementIndex - this->window.size();
  }
  
private:
  input_vector_t& sequence;
  
  AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding;
  
  ComplexHasher hasher;
  
  std::queue<encoded_elem_t> window;
  
  int nextElementIndex;
  
};

#endif
