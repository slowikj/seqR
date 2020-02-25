#ifndef ALPHABET_ENCODER_H
#define ALPHABET_ENCODER_H

#include <functional>
#include <memory>
#include "../dictionary.h"

template<class input_elem_t, class internal_elem_t>
class AlphabetEncoder {
public:
  explicit AlphabetEncoder(std::function<internal_elem_t(input_elem_t)>&& input_to_internal_item_converter)
    : input_to_internal_item_converter(std::move(input_to_internal_item_converter)) { }
  
  template<class encoded_item_t, class input_t>
  Dictionary<internal_elem_t, encoded_item_t> getEncoding(const input_t& input);
  
private:
  std::function<internal_elem_t(input_elem_t)> input_to_internal_item_converter;
  
};

template<class input_elem_t, class internal_elem_t>
template<class encoded_item_t, class input_t>
Dictionary<internal_elem_t, encoded_item_t>
    AlphabetEncoder<input_elem_t, internal_elem_t>::getEncoding(const input_t &input) {
      encoded_item_t currentNum = 1;
      auto res = Dictionary<internal_elem_t, encoded_item_t>();
      for(const input_elem_t& inputElem: input) {
        internal_elem_t internalElem = this->input_to_internal_item_converter(inputElem);
        if(!res.isPresent(internalElem)) {
          res[internalElem] = currentNum++;
        }
      }
      return res;
    }

#endif //ALPHABET_ENCODER_H
