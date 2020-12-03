#ifndef SEQR_ITERATOR_H
#define SEQR_ITERATOR_H

#include <iterator>

namespace dictionary {
    template<class ITER_V, class INNER_CONTAINER_IT_T>
    class iterator_t {
    public:
        // iterator_traits
        using value_type = ITER_V;
        using pointer = ITER_V *;
        using reference = ITER_V &;
        using iterator_category = std::forward_iterator_tag;
        using difference_type = void;

        explicit iterator_t(const INNER_CONTAINER_IT_T &it) : container_iterator_(it) {}

        inline iterator_t<ITER_V, INNER_CONTAINER_IT_T> &operator++() {
            ++this->container_iterator_;
            return *this;
        }

        inline iterator_t<ITER_V, INNER_CONTAINER_IT_T> operator++(int) {
            auto old = iterator(this->container_iterator_);
            ++this->container_iterator_;
            return old;
        }

        inline value_type operator*() {
            return *(this->container_iterator_);
        }

        inline iterator_t<ITER_V, INNER_CONTAINER_IT_T> operator->() {
            return (this->container_iterator_);
        }

        inline bool operator==(const iterator_t<ITER_V, INNER_CONTAINER_IT_T> &other) const {
            return this->container_iterator_ == other.container_iterator_;
        }

        inline bool operator!=(const iterator_t<ITER_V, INNER_CONTAINER_IT_T> &other) const {
            return !((*this) == other);
        }

    private:
        INNER_CONTAINER_IT_T container_iterator_;

    };
}


#endif //SEQR_ITERATOR_H
