#ifndef ROLLING_WINDOW_H
#define ROLLING_WINDOW_H

#include <memory>
#include "../alphabet_encoder/alphabet_encoder.h"
#include "../hash/complex_hasher.h"

template<class input_vector_t, class input_elem_t, class encoded_elem_t>
class RollingWindow {
public:
  
  RollingWindow(input_vector_t& sequence,
                std::unique_ptr<ComplexHasher>&& hasher,
                std::unique_ptr<AlphabetEncoding<input_vector_t, input_elem_t, encoded_elem_t>>&& alphabetEncoding):
    sequence(sequence),
    hasher(std::move(hasher)),
    alphabetEncoding(std::move(alphabetEncoding)) {
    this->nextElementIndex = 0;
  }
  
  void resetIndex(int nextElementIndex) {
    this->nextElementIndex = nextElementIndex;
    this->window.clear();
    this->hasher->clear();
  }
  
  void append() {
    encoded_elem_t encodedElem = this->alphabetEncoding->encode(
      this->sequence[this->nextElementIndex]
    );
    this->window.push(encodedElem);
    this->hasher->append(elem);
    ++this->nextElementIndex;
  }
  
  void moveWindowRight() {
    removeFirst();
    append();
  }
  
  void removeFirst() {
    this->hasher->removeFirst(this->window.front());
  }
  
  std::size_t sequenceSize() const {
    return this->sequence.size();
  }
  
  std::vector<int> getWindowedHashes() const {
    return this->hasher.getHashes();
  }
  
  std::vector<int> getWindowedPositionedHashes() const {
    return this->hasher->getHashes(this->nextElementIndex - 1);
  }
  
private:
  input_vector_t& sequence;
  
  std::unique_ptr<AlphabetEncoding<input_vector_t, input_elem_t, encoded_elem_t>> alphabetEncoding;
  
  std::unique_ptr<ComplexHasher> hasher;
  
  std::queue<encoded_elem_t> window;
  
  int nextElementIndex;
  
};

#endif
