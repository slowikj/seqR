#ifndef ALPHABET_ENCODER_H
#define ALPHABET_ENCODER_H

#include <Rcpp.h>
#include <functional>
#include <memory>

namespace alphabetEncoding {

    template<class input_elem_t_, class encoded_elem_t_,
            template<typename input_t, typename encoded_t, typename...> class dictionary_t>
    class DefaultAlphabetEncoder {
    public:
        using encoded_elem_t = encoded_elem_t_;
        using input_elem_t = input_elem_t_;

        DefaultAlphabetEncoder(dictionary_t<input_elem_t, encoded_elem_t> &&encoder,
                               encoded_elem_t notAllowedEncodingNum) :
                encoder(std::move(encoder)),
                notAllowedEncodingNum(notAllowedEncodingNum) {
        }

        DefaultAlphabetEncoder() = default;

        DefaultAlphabetEncoder(
                DefaultAlphabetEncoder<input_elem_t, encoded_elem_t, dictionary_t> &&other) noexcept = default;

        DefaultAlphabetEncoder(const DefaultAlphabetEncoder<input_elem_t, encoded_elem_t, dictionary_t> &) = delete;

        DefaultAlphabetEncoder &
        operator=(const DefaultAlphabetEncoder<input_elem_t, encoded_elem_t, dictionary_t> &) = delete;

        DefaultAlphabetEncoder &
        operator=(DefaultAlphabetEncoder<input_elem_t, encoded_elem_t, dictionary_t> &&other) noexcept = default;

        inline encoded_elem_t encode(const input_elem_t &inputElem) {
            return isAllowed(inputElem) ?
                   this->encoder[inputElem] :
                   getNotAllowedEncodingNum();
        }

        inline encoded_elem_t encodeUnsafe(const input_elem_t &inputElem) {
            return this->encoder[inputElem];
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
    DefaultAlphabetEncoder<input_elem_t, encoded_elem_t, dictionary_t> getDefaultAlphabetEncoder(input_t &input) {
        encoded_elem_t currentNum = 2;
        dictionary_t<input_elem_t, encoded_elem_t> encoder;
        for (const auto &inputElem: input) {
            if (!encoder.isPresent(inputElem)) {
                encoder[inputElem] = currentNum++;
            }
        }
        return DefaultAlphabetEncoder<input_elem_t, encoded_elem_t, dictionary_t>(
                std::move(encoder),
                1
        );
    }
}

#endif
