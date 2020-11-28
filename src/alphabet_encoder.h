#ifndef ALPHABET_ENCODER_H
#define ALPHABET_ENCODER_H

#include <Rcpp.h>
#include <functional>
#include <memory>

template<class input_elem_t, class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class dictionary_t>
class AlphabetEncoding {
public:
    AlphabetEncoding(dictionary_t<input_elem_t, encoded_elem_t> &&encoder,
                     encoded_elem_t notAllowedEncodingNum) :
            encoder(std::move(encoder)),
            notAllowedEncodingNum(notAllowedEncodingNum) {
    }

    AlphabetEncoding() = default;

    AlphabetEncoding(AlphabetEncoding<input_elem_t, encoded_elem_t, dictionary_t> &&other) noexcept = default;

    AlphabetEncoding(const AlphabetEncoding<input_elem_t, encoded_elem_t, dictionary_t> &) = delete;

    AlphabetEncoding &operator=(const AlphabetEncoding<input_elem_t, encoded_elem_t, dictionary_t> &) = delete;

    AlphabetEncoding &
    operator=(AlphabetEncoding<input_elem_t, encoded_elem_t, dictionary_t> &&other) noexcept = default;

    inline encoded_elem_t encode(const input_elem_t &inputElem) {
        return isAllowed(inputElem) ?
               this->encoder[inputElem] :
               getNotAllowedEncodingNum();
    }

    inline bool isAllowed(const input_elem_t &inputElem) const {
        return this->encoder.isPresent(inputElem);
    }

    inline std::size_t size() const {
        return this->encoder.size();
    }

private:
    dictionary_t<input_elem_t, encoded_elem_t> encoder;
    encoded_elem_t notAllowedEncodingNum;

    inline encoded_elem_t getNotAllowedEncodingNum() const {
        return this->notAllowedEncodingNum;
    }
};

template<class input_t, class input_elem_t, class encoded_elem_t,
        template<typename key, typename value, typename...> class dictionary_t>
inline
AlphabetEncoding<input_elem_t, encoded_elem_t, dictionary_t> getAlphabetEncoding(input_t &input) {
    encoded_elem_t currentNum = 2;
    dictionary_t<input_elem_t, encoded_elem_t> encoder;
    for (const auto &inputElem: input) {
        if (!encoder.isPresent(inputElem)) {
            encoder[inputElem] = currentNum++;
        }
    }
    return AlphabetEncoding<input_elem_t, encoded_elem_t, dictionary_t>(
            std::move(encoder),
            1
    );
}

#endif
