#ifndef SECOND_SINGLE_HASHER_H
#define SECOND_SINGLE_HASHER_H

#include<queue>
#include<functional>
#include<utility>
#include "globals.h"

namespace hashing {

    class SingleHasher {
    public:
        using elem_t = uint32_t;
        using hash_t = config::single_hash_t;

        virtual void append(elem_t elem) = 0;

        virtual void removeFirst(elem_t elem) = 0;

        virtual hash_t getHash() const {
            return currentHash;
        }

        virtual void clear() {
            this->currentHash = 0;
        }

        virtual ~SingleHasher() = default;

    protected:
        hash_t currentHash = 0;

    };
}

#endif
