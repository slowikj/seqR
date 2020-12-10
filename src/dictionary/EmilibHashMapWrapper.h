#ifndef SEQR_EMILIBHASHMAPWRAPPER_H
#define SEQR_EMILIBHASHMAPWRAPPER_H

#include <functional>
#include "dict_iterator.h"
#include "../../inst/include/emilib/hash_map.hpp"

namespace dictionary {
    template<class K, class V, class Hash=std::hash<K>>
    class EmilibHashMapWrapper {
    public:
        using iterator = iterator_t<std::pair<const K &, V>, typename emilib::HashMap<K, V, Hash>::iterator>;
        using const_iterator = iterator_t<const std::pair<const K &, V>, typename emilib::HashMap<K, V, Hash>::const_iterator>;

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

        EmilibHashMapWrapper() = default;

        EmilibHashMapWrapper(const EmilibHashMapWrapper &) = default;

        EmilibHashMapWrapper &operator=(const EmilibHashMapWrapper &) = default;

        EmilibHashMapWrapper(EmilibHashMapWrapper &&) noexcept = default;

        EmilibHashMapWrapper &operator=(EmilibHashMapWrapper &&) noexcept = default;

        ~EmilibHashMapWrapper() = default;

    private:
        emilib::HashMap<K, V, Hash> inner_map_;
    };
}

#endif //SEQR_EMILIBHASHMAPWRAPPER_H
