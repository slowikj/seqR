#ifndef SEQR_ORDERED_MAP_WRAPPER_H
#define SEQR_ORDERED_MAP_WRAPPER_H

#include <map>

namespace dictionary {

    template<class K, class V, class...>
    class OrderedMapWrapper {
    public:
        using iterator = typename std::map<K, V>::iterator;
        using const_iterator = typename std::map<K, V>::const_iterator;

        inline V &operator[](const K &key) {
            return this->inner_map_[key];
        }

        inline V &operator[](K &&key) {
            return this->inner_map_[std::move(key)];
        }

        inline std::size_t size() const {
            return this->inner_map_.size();
        }

        inline bool isPresent(const K &key) const {
            return this->inner_map_.find(key) != std::end(this->inner_map_);
        }

        inline iterator begin() {
            return iterator(this->inner_map_.begin());
        }

        inline const_iterator begin() const noexcept {
            return const_iterator(this->inner_map_.begin());
        }

        inline iterator end() {
            return iterator(this->inner_map_.end());
        }

        inline const_iterator end() const noexcept {
            return const_iterator(this->inner_map_.end());
        }

        OrderedMapWrapper() = default;

        OrderedMapWrapper(const OrderedMapWrapper &) = default;

        OrderedMapWrapper &operator=(const OrderedMapWrapper &) = default;

        OrderedMapWrapper(OrderedMapWrapper &&) noexcept = default;

        OrderedMapWrapper &operator=(OrderedMapWrapper &&) noexcept = default;

        ~OrderedMapWrapper() = default;

    private:
        std::map<K, V> inner_map_;
    };
}

#endif //SEQR_ORDERED_MAP_WRAPPER_H
