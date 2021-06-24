#pragma once

#include <vector>

#include "globals.h"
#include "polynomial_single_hasher.h"

namespace hashing {

template <class encoded_sequence_t>
class PrefixSequencePolynomialHasher {
 public:
  using hash_t = config::multidim_hash_t;

  PrefixSequencePolynomialHasher(const encoded_sequence_t &sequence,
                                 const std::vector<PolynomialSingleHasherConfig> &polynomialHasherConfigs)
      : polynomialHasherConfigs(polynomialHasherConfigs) {
    initModuloMComputers(polynomialHasherConfigs);
    computePrefixValues(sequence);
  }

  [[nodiscard]] inline hash_t getHash(std::size_t begin, std::size_t end) const {
    hash_t res(polynomialHasherConfigs.size());
    for (std::size_t hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
      auto M = polynomialHasherConfigs[hasherInd].M;
      res[hasherInd] = moduloMComputers[hasherInd].get(
          M + prefixComplexHashes[end + 1][hasherInd] -
          moduloMComputers[hasherInd].get(
              prefixComplexHashes[begin][hasherInd] * prefixP[end - begin + 1][hasherInd]));
    }
    return res;
  }

  [[nodiscard]] inline hash_t getHashForSeveralIntervals(
      std::size_t beginPosition,
      const std::vector<std::pair<std::size_t, std::size_t>> &contiguousIntervals) const {
    hash_t res(this->getHashersNum());
    for (const auto &interval : contiguousIntervals) {
      auto intervalHash = this->getHash(
          interval.first + beginPosition,
          interval.second + beginPosition);
      std::size_t intervalLength = util::getIntervalLength(interval);
      for (std::size_t hasherInd = 0; hasherInd < res.size(); ++hasherInd) {
        int powerP = this->getHasherP(hasherInd, intervalLength);
        res[hasherInd] = moduloMComputers[hasherInd].get(
            res[hasherInd] * powerP + intervalHash[hasherInd]);
      }
    }
    return res;
  }

 private:
  const std::vector<PolynomialSingleHasherConfig> &polynomialHasherConfigs;
  std::vector<hash_t> prefixP;
  std::vector<hash_t> prefixComplexHashes;
  std::vector<util::ModuloComputer> moduloMComputers;

  void initModuloMComputers(const std::vector<PolynomialSingleHasherConfig> &configs) {
    for (const auto &config : configs) {
      moduloMComputers.emplace_back(config.M);
    }
  }

  inline void computePrefixValues(const encoded_sequence_t &sequence) {
    initPrefixP(sequence.size(), polynomialHasherConfigs.size());
    initPrefixComplexHashes(sequence.size(), polynomialHasherConfigs.size());
    for (std::size_t seqInd = 0; seqInd < sequence.size(); ++seqInd) {
      appendPrefixValues(sequence, seqInd);
    }
  }

  inline void initPrefixP(std::size_t sequenceLength, std::size_t hashersNum) {
    prefixP.reserve(sequenceLength);
    prefixP.push_back(hash_t(hashersNum, 1));
  }

  inline void initPrefixComplexHashes(std::size_t sequenceLength, std::size_t hashersNum) {
    prefixComplexHashes.reserve(sequenceLength);
    prefixComplexHashes.push_back(hash_t(hashersNum));
  }

  inline void appendPrefixValues(const encoded_sequence_t &sequence,
                                 std::size_t seqInd) {
    appendCurrentComplexHash(sequence[seqInd]);
    appendCurrentPowerP();
  }

  inline void appendCurrentComplexHash(const typename encoded_sequence_t::encoded_elem_t &encodedElem) {
    hash_t prefixHash(polynomialHasherConfigs.size());
    for (std::size_t hasherInd = 0; hasherInd < prefixHash.size(); ++hasherInd) {
      auto P = polynomialHasherConfigs[hasherInd].P;
      prefixHash[hasherInd] = moduloMComputers[hasherInd].get(
          prefixComplexHashes.back()[hasherInd] * P + encodedElem);
    }
    prefixComplexHashes.push_back(std::move(prefixHash));
  }

  inline void appendCurrentPowerP() {
    hash_t powersP(polynomialHasherConfigs.size());
    for (std::size_t hasherInd = 0; hasherInd < powersP.size(); ++hasherInd) {
      auto P = polynomialHasherConfigs[hasherInd].P;
      powersP[hasherInd] = moduloMComputers[hasherInd].get(prefixP.back()[hasherInd] * P);
    }
    prefixP.push_back(std::move(powersP));
  }

  [[nodiscard]] inline PolynomialSingleHasherConfig::elem_t getHasherP(std::size_t hasherIndex, std::size_t power = 1) const {
    return this->prefixP[power][hasherIndex];
  }

  [[nodiscard]] inline PolynomialSingleHasherConfig::elem_t getHasherM(std::size_t hasherIndex) const {
    return this->polynomialHasherConfigs[hasherIndex].M;
  }

  [[nodiscard]] inline std::size_t getHashersNum() const {
    return this->polynomialHasherConfigs.size();
  }
};
}  // namespace hashing
