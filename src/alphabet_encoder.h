#ifndef FIRST_ALPHABET_ENCODER_H
#define FIRST_ALPHABET_ENCODER_H

#include<functional>
#include "dictionary/dictionary.h"
#include "dictionary/unordered_map_dictionary.h"

template<class input_elem_t, class internal_elem_t>
class alphabet_encoder {
public:
  template<class encoded_item_t, class input_t>
  dictionary<internal_elem_t, encoded_item_t> get_encoder(const input_t& input);
  
private:
  const std::function<internal_elem_t(input_elem_t)> &input2internal_item_converter;
  
};

template<class input_elem_t, class internal_elem_t>
template<class encoded_item_t, class input_t>
dictionary<internal_elem_t, encoded_item_t>
alphabet_encoder<input_elem_t, internal_elem_t>::get_encoder(const input_t &input) {
  int current_num = 1;
  auto res = unordered_map_dictionary<internal_elem_t, encoded_item_t>();
  for(const input_elem_t& inputElem: input) {
    internal_elem_t internalElem = this->input2internal_item_converter(inputElem);
    if(!res.is_present(internalElem)) {
      res[internalElem] = current_num++;
    }
  }
  return res;
}


#endif //FIRST_ALPHABET_ENCODER_H
