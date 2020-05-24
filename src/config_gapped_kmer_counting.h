#ifndef SEQR_CONFIG_GAPPED_KMER_COUNTING_H
#define SEQR_CONFIG_GAPPED_KMER_COUNTING_H

#include <vector>
#include "hash/polynomial_single_hasher.h"

inline std::vector<PolynomialSingleHasherConfig> getGappedKMerHasherConfigs() {
    std::vector<PolynomialSingleHasherConfig> res;
    res.emplace_back(101, 1e9 + 7);
    res.emplace_back(97, 1e9 + 33);
    return std::move(res);
}

#endif //SEQR_CONFIG_GAPPED_KMER_COUNTING_H