#pragma once

#include <utility>
#include <memory>

template <class encoded_elem_t>
class EncodedSequenceRow {
public:
    EncodedSequenceRow(std::shared_ptr<encoded_elem_t[]> encoded, int begin, int size)
            : encoded_(std::move(encoded)), begin_(begin), size_(size) {}

    inline const encoded_elem_t &operator[](int index) const {
        return this->encoded_[begin_ + index];
    }

    inline std::size_t size() const {
        return this->size_;
    }

private:
    std::shared_ptr<encoded_elem_t[]> encoded_;
    std::size_t begin_, size_;
};