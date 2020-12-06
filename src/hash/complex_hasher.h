#ifndef SECOND_COMPLEX_HASHER_H
#define SECOND_COMPLEX_HASHER_H

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>
#include "single_hasher.h"
#include "types.h"

namespace hashing {

    class ComplexHasher {
    public:
        explicit ComplexHasher(std::vector<std::unique_ptr<SingleHasher>> &&singleHashers) :
                singleHashers(std::move(singleHashers)) {
        }

        inline void clear() {
            std::for_each(std::begin(this->singleHashers), std::end(this->singleHashers),
                          [](std::unique_ptr<SingleHasher> &singleHasher) {
                              singleHasher->clear();
                          }
            );
        }

        inline void append(const int &elem) {
            std::for_each(std::begin(this->singleHashers), std::end(this->singleHashers),
                          [&elem](std::unique_ptr<SingleHasher> &singleHasher) {
                              singleHasher->append(elem);
                          });
        }

        inline void removeFirst(const int &elem) {
            std::for_each(std::begin(this->singleHashers), std::end(this->singleHashers),
                          [&elem](std::unique_ptr<SingleHasher> &singleHasher) {
                              singleHasher->removeFirst(elem);
                          });
        }

        [[nodiscard]] inline std::vector<single_hash_t> getHashes(int position) const {
            auto res = getHashes();
            res.push_back(position);
            return res;
        }

        [[nodiscard]] inline std::vector<single_hash_t> getHashes() const {
            return prepareResultHashes(
                    [](const std::unique_ptr<SingleHasher> &singleHasher) -> int {
                        return singleHasher->getHash();
                    }
            );
        }

    private:
        std::vector<std::unique_ptr<SingleHasher>> singleHashers;

        inline std::vector<single_hash_t> prepareResultHashes(
                std::function<single_hash_t(const std::unique_ptr<SingleHasher> &)> &&transformFunc) const {
            std::vector<single_hash_t> resultHashes;
            std::transform(std::begin(singleHashers), std::end(singleHashers),
                           std::back_inserter(resultHashes),
                           transformFunc);
            return resultHashes;
        }

    };
}

#endif
