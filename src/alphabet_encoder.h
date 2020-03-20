#ifndef ALPHABET_ENCODER_H
#define ALPHABET_ENCODER_H

// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>

#include <functional>
#include <memory>
#include "dictionary.h"

template<class input_elem_t, class internal_elem_t, class encoded_elem_t>
class AlphabetEncoding {
public:
  AlphabetEncoding(std::function<internal_elem_t(const input_elem_t&)> inputToInternalItemConverter,
                   Dictionary<internal_elem_t, encoded_elem_t>&& internalToEncoded,
                   encoded_elem_t notAllowedEncodingNum):
    internalToEncoded(std::move(internalToEncoded)),
    inputToInternalItemConverter(inputToInternalItemConverter),
    notAllowedEncodingNum(notAllowedEncodingNum) {
  }
  
  AlphabetEncoding() = default;
  
  AlphabetEncoding(AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>&& other) noexcept = default;
  
  AlphabetEncoding(const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>&) = delete;
  
  AlphabetEncoding& operator=(const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>&) = delete;
  
  AlphabetEncoding& operator=(AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>&& other) noexcept = default;
  
  encoded_elem_t encode(const input_elem_t& inputElem) {
    return internalToEncoded[inputToInternalItemConverter(inputElem)];
  }
  
  bool isAllowed(const input_elem_t& inputElem) const {
    return internalToEncoded.isPresent(inputToInternalItemConverter(inputElem));
  }
  
  encoded_elem_t getNotAllowedEncodingNum() const {
    return this->notAllowedEncodingNum;
  }
  
  std::size_t alphabetSize() const {
    return internalToEncoded.size();
  }
  
private:
  Dictionary<internal_elem_t, encoded_elem_t> internalToEncoded;
  std::function<internal_elem_t(const input_elem_t&)> inputToInternalItemConverter;
  encoded_elem_t notAllowedEncodingNum;
};

template<class input_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t> getAlphabetEncoding(
    input_t& input,
    std::function<internal_elem_t(const input_elem_t&)> inputToInternalItemConverter) {
  encoded_elem_t currentNum = 2;
  Dictionary<internal_elem_t, encoded_elem_t> internalToEncoded;
  for(const auto& inputElem: input) {
    internal_elem_t internalElem = inputToInternalItemConverter(inputElem);
    if(!internalToEncoded.isPresent(internalElem)) {
      internalToEncoded[internalElem] = currentNum++;
    }
  }
  return AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>(
      inputToInternalItemConverter,
      std::move(internalToEncoded),
      1);
}

#endif //ALPHABET_ENCODER_H
