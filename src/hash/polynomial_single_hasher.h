#ifndef SECOND_POLYNOMIAL_SINGLE_HASHER_H
#define SECOND_POLYNOMIAL_SINGLE_HASHER_H

#include "single_hasher.h"
#include "../utils.h"

struct PolynomialSingleHasherConfig {
  
  int P;
  int M;
  
  PolynomialSingleHasherConfig(int P, int M):
    P(P), M(M) {
  }
  
  PolynomialSingleHasherConfig() = delete;
  
  PolynomialSingleHasherConfig(const PolynomialSingleHasherConfig&) = default;
  
  PolynomialSingleHasherConfig& operator=(const PolynomialSingleHasherConfig&) = default;
  
};

class PolynomialSingleHasher : public SingleHasher {
public:
  PolynomialSingleHasher(PolynomialSingleHasherConfig&& config) :
    config(std::move(config)) {
    this->P_M_2 = compute_power_fast(config.P, config.M - 2, config.M);
    this->nextPowerP = 1;
    this->currentPowerP = 0;
  }
  
  void append(const int &elem) override {
    this->currentHash = this->computeHash(this->currentHash, elem);
    this->currentPowerP = this->nextPowerP;
    this->nextPowerP = this->computeNextPowerP(this->nextPowerP);
  }
  
  void removeFirst(const int &elem) override {
    this->currentHash = static_cast<int>(
      (this->currentHash -
        (static_cast<long long>(elem) * this->computePreviousPowerP(this->nextPowerP)) + config.M) % config.M
    );
    this->nextPowerP = this->computePreviousPowerP(nextPowerP);
    this->currentPowerP = this->computePreviousPowerP(nextPowerP);
  }
  
  [[nodiscard]] int getHash() const override {
    return SingleHasher::getHash();
  }
  
  [[nodiscard]] int getHash(int position) const override {
    return this->computeHash(currentHash, position);
  }
  
  void clear() override {
    SingleHasher::clear();
  }
  
  int getP() const {
    return this->config.P;
  }
  
  int getM() const {
    return this->config.M;
  }
  
  int getCurrentPowerP() const {
    return this->currentPowerP;
  }
  
private:
  PolynomialSingleHasherConfig config;
  int P_M_2; // P^(M-2) MOD M
  int nextPowerP;
  int currentPowerP;
  
  int computeHash(int currentHash, const int &elem) const {
    return static_cast<int>(
      (static_cast<long long>(this->currentHash) * config.P + elem) % config.M
    );
  }
  
  int computeNextPowerP(int currentPowerP) const {
    return static_cast<int>(
      (static_cast<long long>(currentPowerP) * config.P) % config.M
    );
  }
  
  int computePreviousPowerP(int currentPowerP) const {
    return static_cast<int>(
      (static_cast<long long>(currentPowerP) * this->P_M_2) % config.M
    );
  }
};

#endif //SECOND_POLYNOMIAL_SINGLE_HASHER_H
