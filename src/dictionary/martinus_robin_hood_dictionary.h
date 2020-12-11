#ifndef SEQR_MARTINUS_ROBIN_HOOD_DICTIONARY_H
#define SEQR_MARTINUS_ROBIN_HOOD_DICTIONARY_H

#include "../../inst/thirdparty/robin_hood_martinus.h"
#include <map>

namespace dictionary {

    template<class K, class V, class Hash=std::hash<K>>
    class MartinusRobinHoodDictionary {
    public:
        using iterator = typename robin_hood::unordered_map<K, V, Hash>::iterator;
        using const_iterator = typename robin_hood::unordered_map<K, V, Hash>::const_iterator ;

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

        MartinusRobinHoodDictionary() = default;

        MartinusRobinHoodDictionary(const MartinusRobinHoodDictionary &) = default;

        MartinusRobinHoodDictionary &operator=(const MartinusRobinHoodDictionary &) = default;

        MartinusRobinHoodDictionary(MartinusRobinHoodDictionary &&) noexcept = default;

        MartinusRobinHoodDictionary &operator=(MartinusRobinHoodDictionary &&) noexcept = default;

        ~MartinusRobinHoodDictionary() = default;

    private:
        robin_hood::unordered_map<K, V, Hash> inner_map_;
    };
}

#endif //SEQR_MARTINUS_ROBIN_HOOD_DICTIONARY_H
