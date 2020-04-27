#ifndef ALPHABET_ENCODER_H
#define ALPHABET_ENCODER_H

#include <Rcpp.h>
#include <functional>
#include <memory>
#include "dictionary.h"

template<class input_elem_t, class encoded_elem_t, class hasher_t>
class AlphabetEncoding {
public:
    AlphabetEncoding(Dictionary<input_elem_t, encoded_elem_t, hasher_t> &&encoder,
                     encoded_elem_t notAllowedEncodingNum) :
            encoder(std::move(encoder)),
            notAllowedEncodingNum(notAllowedEncodingNum) {
    }

    AlphabetEncoding() = default;

    AlphabetEncoding(AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t> &&other) noexcept = default;

    AlphabetEncoding(const AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t> &) = delete;

    AlphabetEncoding &operator=(const AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t> &) = delete;

    AlphabetEncoding &operator=(AlphabetEncoding<input_elem_t, encoded_elem_t, hasher_t> &&other) noexcept = default;

    inline encoded_elem_t encode(const input_elem_t &inputElem) {
        return this->encoder[inputElem];
    }

    inline bool isAllowed(const input_elem_t &inputElem) const {
        return this->encoder.isPresent(inputElem);
    }

    inline encoded_elem_t getNotAllowedEncodingNum() const {
        return this->notAllowedEncodingNum;
    }

    inline std::size_t size() const {
        return this->encoder.size();
    }

private:
    Dictionary<input_elem_t, encoded_elem_t, hasher_t> encoder;
    encoded_elem_t notAllowedEncodingNum;
};

template<class input_t, class input_elem_t, class encoded_elem_t, class alphabet_hasher_t>
inline
AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> getAlphabetEncoding(input_t &input) {
    encoded_elem_t currentNum = 2;
    Dictionary<input_elem_t, encoded_elem_t, alphabet_hasher_t> encoder;
    for (const auto &inputElem: input) {
        if (!encoder.isPresent(inputElem)) {
            encoder[inputElem] = currentNum++;
        }
    }
    return AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t>(
            std::move(encoder),
            1
    );
}

template<class encoded_elem_t, class alphabet_hasher_t>
inline
AlphabetEncoding<encoded_elem_t, encoded_elem_t, alphabet_hasher_t> prepareAlphabetEncodingForTidysq(
        Rcpp::StringVector &alphabet,
        Rcpp::StringVector &elementsEncoding) {
    Dictionary<encoded_elem_t, encoded_elem_t, alphabet_hasher_t> encoder;
    for (const auto &alphabetElem: alphabet) {
        for (int encoding_i = 0; encoding_i < elementsEncoding.size(); ++encoding_i) {
            if (alphabetElem == elementsEncoding[encoding_i]) {
                auto index = static_cast<encoded_elem_t>(encoding_i + 1);
                encoder[index] = index;
            }
        }
    }
    return AlphabetEncoding<encoded_elem_t, encoded_elem_t, alphabet_hasher_t>(
            std::move(encoder),
            elementsEncoding.size() + 1
    );
}

#endif
