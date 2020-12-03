#ifndef SEQR_PRIMES_H
#define SEQR_PRIMES_H

#include <vector>
#include <utility>

namespace hashing {

    const static std::vector<std::pair<int, int>> hashPrimes = {
            {101, 1e9 + 7},
            {97, 1e9 + 33},
            {89, 1e9 + 9},
            {71, 1192901131},
            {67, 1190705431},
            {61, 1190494969},
            {59, 1086218491},
            {53, 1065433423}
    };
}

#endif //SEQR_PRIMES_H
