#pragma once

namespace dictionary {

    template<class K, class V,
            template <typename K_, typename V_, typename...> class container_t,
            class ...WArgs>
    class StlLikeDictionaryWrapper {
    public:
        using iterator = typename container_t<K, V, WArgs...>::iterator;
        using const_iterator = typename container_t<K, V, WArgs...>::const_iterator;

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

        StlLikeDictionaryWrapper() = default;

        StlLikeDictionaryWrapper(const StlLikeDictionaryWrapper &) = default;

        StlLikeDictionaryWrapper &operator=(const StlLikeDictionaryWrapper &) = default;

        StlLikeDictionaryWrapper(StlLikeDictionaryWrapper
        &&) = default;

        StlLikeDictionaryWrapper &operator=(StlLikeDictionaryWrapper &&) = default;

        ~StlLikeDictionaryWrapper() = default;

    private:
        container_t<K, V, WArgs...> inner_map_;
    };

}
