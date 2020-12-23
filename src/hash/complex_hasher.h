#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>
#include "single_hasher.h"
#include "globals.h"

namespace hashing {

    class ComplexHasher {
    public:
        using elem_t = SingleHasher::elem_t;
        using hash_t = config::multidim_hash_t;

        explicit ComplexHasher(std::vector<std::unique_ptr<SingleHasher>> &&singleHashers) :
                singleHashers(std::move(singleHashers)) {
        }

        inline void clear() {
            std::for_each(
                    std::begin(this->singleHashers),
                    std::end(this->singleHashers),
                    [](std::unique_ptr<SingleHasher> &singleHasher) {
                        singleHasher->clear();
                    }
            );
        }

        inline void append(elem_t elem) {
            std::for_each(
                    std::begin(this->singleHashers),
                    std::end(this->singleHashers),
                    [&elem](std::unique_ptr<SingleHasher> &singleHasher) {
                        singleHasher->append(elem);
                    });
        }

        inline void removeFirst(elem_t elem) {
            std::for_each(
                    std::begin(this->singleHashers),
                    std::end(this->singleHashers),
                    [&elem](std::unique_ptr<SingleHasher> &singleHasher) {
                        singleHasher->removeFirst(elem);
                    });
        }

        [[nodiscard]] inline hash_t getHashes(std::size_t position) const {
            auto res = getHashes();
            res.push_back(position);
            return res;
        }

        [[nodiscard]] inline hash_t getHashes() const {
            return prepareResultHashes(
                    [](const std::unique_ptr<SingleHasher> &singleHasher) -> config::single_hash_t {
                        return singleHasher->getHash();
                    }
            );
        }

    private:
        std::vector<std::unique_ptr<SingleHasher>> singleHashers;

        inline hash_t prepareResultHashes(
                std::function<config::single_hash_t(const std::unique_ptr<SingleHasher> &)> &&transformFunc) const {
            hash_t resultHashes;
            std::transform(
                    std::begin(singleHashers),
                    std::end(singleHashers),
                    std::back_inserter(resultHashes),
                    transformFunc);
            return resultHashes;
        }
    };
}
