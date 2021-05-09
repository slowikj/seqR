#pragma once

#include <iterator>

namespace dictionary {
template <class iter_value_t, class inner_container_iterator_t>
class iterator_t {
 public:
  // iterator_traits
  using value_type = iter_value_t;
  using pointer = iter_value_t *;
  using reference = iter_value_t &;
  using iterator_category = std::forward_iterator_tag;
  using difference_type = void;

  explicit iterator_t(const inner_container_iterator_t &it) : container_iterator_(it) {}

  inline iterator_t<iter_value_t, inner_container_iterator_t> &operator++() {
    ++this->container_iterator_;
    return *this;
  }

  inline iterator_t<iter_value_t, inner_container_iterator_t> operator++(int) {
    auto old = iterator(this->container_iterator_);
    ++this->container_iterator_;
    return old;
  }

  inline value_type operator*() {
    return *(this->container_iterator_);
  }

  inline iterator_t<iter_value_t, inner_container_iterator_t> operator->() {
    return (this->container_iterator_);
  }

  inline bool operator==(const iterator_t<iter_value_t, inner_container_iterator_t> &other) const {
    return this->container_iterator_ == other.container_iterator_;
  }

  inline bool operator!=(const iterator_t<iter_value_t, inner_container_iterator_t> &other) const {
    return !((*this) == other);
  }

 private:
  inner_container_iterator_t container_iterator_;
};
}  // namespace dictionary
