#ifndef SEQR_LINEAR_LIST_DICTIONARY_H
#define SEQR_LINEAR_LIST_DICTIONARY_H

#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include "dict_iterator.h"

namespace dictionary {
    template<class K, class V, class...>
    class LinearListDictionary {
    public:
        using iterator = iterator_t<std::pair<K, V>, typename std::vector<std::pair<K, V>>::iterator>;
        using const_iterator = iterator_t<const std::pair<K, V>, typename std::vector<std::pair<K, V>>::const_iterator>;

        LinearListDictionary() = default;

        LinearListDictionary(const LinearListDictionary &) = default;

        LinearListDictionary &operator=(const LinearListDictionary &) = default;

        LinearListDictionary &operator=(LinearListDictionary &&) noexcept = default;

        LinearListDictionary(LinearListDictionary &&) noexcept = default;

        inline V &operator[](const K &key) {
            return getValueOrNewlyAddedDefault(key);
        }

        inline V &operator[](K &&key) {
            return getValueOrNewlyAddedDefault(key);
        }

        inline std::size_t size() const {
            return items.size();
        }

        inline bool isPresent(const K &key) {
            return std::any_of(
                    std::begin(items),
                    std::end(items),
                    [this, &k = std::as_const(key)](const std::pair<K, V> &item) {
                        return item.first == k;
                    }
            );
        }

        iterator begin() {
            return iterator(items.begin());
        }

        iterator end() {
            return iterator(items.end());
        }

        const_iterator begin() const {
            return const_iterator(items.begin());
        }

        const_iterator end() const {
            return const_iterator(items.end());
        }

    private:
        std::vector<std::pair<K, V>> items;

        inline V &getValueOrNewlyAddedDefault(const K &key) {
            for (auto &pair: items) {
                if (pair.first == key) {
                    return pair.second;
                }
            }
            items.emplace_back(key, std::move(V()));
            return items[items.size() - 1].second;
        }

    };

}
#endif //SEQR_LINEAR_LIST_DICTIONARY_H
