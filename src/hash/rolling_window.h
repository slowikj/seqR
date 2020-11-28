#ifndef ROLLING_WINDOW_H
#define ROLLING_WINDOW_H

#include <memory>
#include "../default_alphabet_encoder.h"
#include "../hash/complex_hasher.h"
#include "../utils.h"

template<class input_vector_t, class alphabet_encoding_t>
class RollingWindow {
public:
    using encoded_elem_t = typename alphabet_encoding_t::encoded_elem_t;

    RollingWindow(input_vector_t &sequence,
                  ComplexHasher &&hasher,
                  alphabet_encoding_t &alphabetEncoding) :
            sequence(sequence),
            hasher(std::move(hasher)),
            alphabetEncoding(alphabetEncoding) {
        this->nextElementIndex = 0;
    }

    RollingWindow() = delete;

    inline void resetIndex(int nextElementIndex) {
        this->nextElementIndex = nextElementIndex;
        clear<encoded_elem_t>(this->window);
        this->hasher.clear();
    }

    inline void append() {
        encoded_elem_t encodedElem = this->alphabetEncoding.encodeUnsafe(
                this->sequence[this->nextElementIndex]
        );
        this->window.push(encodedElem);
        this->hasher.append(encodedElem);
        ++this->nextElementIndex;
    }

    inline void moveWindowRight() {
        removeFirst();
        append();
    }

    inline void removeFirst() {
        this->hasher.removeFirst(this->window.front());
        this->window.pop();
    }

    inline std::size_t sequenceSize() const {
        return this->sequence.size();
    }

    inline std::vector<int> getWindowedHashes() const {
        return this->hasher.getHashes();
    }

    inline std::vector<int> getWindowedPositionedHashes() const {
        return this->hasher.getHashes(this->nextElementIndex);
    }

    inline int currentBeginIndex() const {
        return this->nextElementIndex - this->window.size();
    }

private:
    input_vector_t &sequence;

    alphabet_encoding_t &alphabetEncoding;

    ComplexHasher hasher;

    std::queue<encoded_elem_t> window;

    int nextElementIndex;
};

#endif
