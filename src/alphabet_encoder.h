#ifndef ALPHABET_ENCODER_H
#define ALPHABET_ENCODER_H

// [[Rcpp::plugins("c++17")]]
#include <Rcpp.h>
#include <functional>
#include <memory>
#include "dictionary.h"

template<class input_elem_t, class encoded_elem_t, class hasher_t>
class AlphabetEncoding {
public:
  AlphabetEncoding(Dictionary<input_elem_t, encoded_elem_t, hasher_t>&& encoder,
                   encoded_elem_t notAllowedEncodingNum):
    encoder(std::move(encoder)),
    notAllowedEncodingNum(notAllowedEncodingNum) {
  }
  
  AlphabetEncoding() = default;
  
  AlphabetEncoding(AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t>&& other) noexcept = default;
  
  AlphabetEncoding(const AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t>&) = delete;
  
  AlphabetEncoding& operator=(const AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t>&) = delete;
  
  AlphabetEncoding& operator=(AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t>&& other) noexcept = default;
  
  encoded_elem_t encode(const input_elem_t& inputElem) {
    return this->encoder[inputElem];
  }
  
  bool isAllowed(const input_elem_t& inputElem) const {
    return this->encoder.isPresent(inputElem);
  }
  
  encoded_elem_t getNotAllowedEncodingNum() const {
    return this->notAllowedEncodingNum;
  }
  
  std::size_t alphabetSize() const {
    return this->encoder.size();
  }
  
private:
  Dictionary<input_elem_t, encoded_elem_t, hasher_t> encoder;
  encoded_elem_t notAllowedEncodingNum;
};

template<class input_t, class input_elem_t, class encoded_elem_t, class hasher_t>
AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t> getAlphabetEncoding(input_t& input) {
  encoded_elem_t currentNum = 2;
  Dictionary<input_elem_t, encoded_elem_t, hasher_t> encoder;
  for(const auto& inputElem: input) {
    if(!encoder.isPresent(inputElem)) {
      encoder[inputElem] = currentNum++;
    }
  }
  return AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t>(
      std::move(encoder),
      1
  );
}

#endif
