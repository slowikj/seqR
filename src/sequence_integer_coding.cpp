
// [[Rcpp::plugins("c++1z")]]
#include<Rcpp.h>
#include<string>
#include<unordered_map>
#include<vector>
#include<functional>
#include<tuple>

template<class ALPHABET_VECTOR_TYPE, class CPP_ITEM_TYPE>
std::tuple<std::unordered_map<CPP_ITEM_TYPE, short>*, std::unordered_map<short, CPP_ITEM_TYPE>*>
enumerate_alphabet(ALPHABET_VECTOR_TYPE& alphabet) {
  short least_available_number = 1;
  auto val2short_encoder = new std::unordered_map<CPP_ITEM_TYPE, short>();
  auto short2val_decoder = new std::unordered_map<short, CPP_ITEM_TYPE>();
  for(const auto& alphabet_item: alphabet) {
    CPP_ITEM_TYPE str = Rcpp::as<CPP_ITEM_TYPE>(alphabet_item);
    if(val2short_encoder->find(str) == val2short_encoder->end()) {
      (*val2short_encoder)[str] = least_available_number;
      (*short2val_decoder)[least_available_number] = str;
      ++least_available_number;
    }
  }
  return {val2short_encoder, short2val_decoder};
}

const short NOT_ALLOWED_CHARACTER_CODE = -1;

template<class VECTOR_TYPE, class ITEM_TYPE>
std::vector<short> enumerate_sequence(VECTOR_TYPE& seq,
                                      std::unordered_map<ITEM_TYPE, short>* val2short_encoder) {
  std::vector<short> res;
  res.reserve(seq.size());
  for(const auto& seq_item: seq) {
    ITEM_TYPE item_cpp = Rcpp::as<ITEM_TYPE>(seq_item);
    short item_encoded = val2short_encoder -> find(item_cpp) == val2short_encoder -> end()
      ? NOT_ALLOWED_CHARACTER_CODE
      : (*val2short_encoder)[item_cpp];
    res.push_back(item_encoded);
  }
  return res;
}

template<class VECTOR_TYPE, class ITEM_TYPE>
std::vector<short> enumerate_sequence(VECTOR_TYPE& v,
                                      VECTOR_TYPE& alphabet) {
  auto [val2short_encoder, short2val_decoder] = enumerate_alphabet<VECTOR_TYPE, ITEM_TYPE>(alphabet);
  auto res =enumerate_sequence<VECTOR_TYPE, ITEM_TYPE>(v, val2short_encoder);
  delete val2short_encoder;
  delete short2val_decoder;
  return res;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector enumerate_string_sequence(Rcpp::StringVector& sequence,
                                              Rcpp::StringVector& alphabet) {
  return Rcpp::wrap(enumerate_sequence<Rcpp::StringVector, std::string>(sequence, alphabet));
}
