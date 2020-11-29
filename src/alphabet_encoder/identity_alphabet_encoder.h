#ifndef SEQR_IDENTITY_ALPHABET_ENCODER_H
#define SEQR_IDENTITY_ALPHABET_ENCODER_H

#include <utility>

template<class encoded_elem_t_>
class IdentityAlphabetEncoder {
public:
    using encoded_elem_t = encoded_elem_t_;
    using input_elem_t = encoded_elem_t;

    IdentityAlphabetEncoder(encoded_elem_t_ allowedBegin,
                            encoded_elem_t_ allowedEnd)
            : allowedBegin(allowedBegin), allowedEnd(allowedEnd) {}

    inline encoded_elem_t encode(const input_elem_t &inputElem) {
        return inputElem;
    }

    inline encoded_elem_t encodeUnsafe(const input_elem_t &inputElem) {
        return encode(inputElem);
    }

    inline bool isAllowed(const input_elem_t &inputElem) const {
        return allowedBegin <= inputElem && inputElem <= allowedEnd;
    }

    inline std::size_t size() const {
        return (allowedEnd - allowedBegin + 1);
    }

private:
    encoded_elem_t_ allowedBegin, allowedEnd;

};

#endif //SEQR_IDENTITY_ALPHABET_ENCODER_H
