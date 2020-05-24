#ifndef SEQR_CONFIG_KMER_COUNTING_H
#define SEQR_CONFIG_KMER_COUNTING_H

#include "hash/complex_hasher.h"
#include "hash/polynomial_single_hasher.h"

inline ComplexHasher createKMerComplexHasher() {
    std::vector<std::unique_ptr<SingleHasher>> singleHashers;
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(101, 1e9 + 33)));
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(97, 1e9 + 7)));
    ComplexHasher complexHasher(std::move(singleHashers));
    return complexHasher;
}

#endif //SEQR_CONFIG_KMER_COUNTING_H