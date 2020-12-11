#ifndef SEQR_STL_DICTIONARY_TEMPLATE_H
#define SEQR_STL_DICTIONARY_TEMPLATE_H

namespace dictionary {

    template<class K, class V,
            template <typename K_, typename V_, typename H > class container_t,
            class ...WArgs>
    class StlDictionaryWrapper {
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

        StlDictionaryWrapper() = default;

        StlDictionaryWrapper(const StlDictionaryWrapper &) = default;

        StlDictionaryWrapper &operator=(const StlDictionaryWrapper &) = default;

        StlDictionaryWrapper(StlDictionaryWrapper
        &&) noexcept = default;

        StlDictionaryWrapper &operator=(StlDictionaryWrapper &&) noexcept = default;

        ~StlDictionaryWrapper() = default;

    private:
        container_t<K, V, WArgs...> inner_map_;
    };

}

#endif //SEQR_STL_DICTIONARY_TEMPLATE_H
