#ifndef SECOND_POLYNOMIAL_SINGLE_HASHER_H
#define SECOND_POLYNOMIAL_SINGLE_HASHER_H

#include "single_hasher.h"
#include "../utils.h"

class PolynomialSingleHasher : public SingleHasher {
public:
  PolynomialSingleHasher(int P, int M) :
    P(P), M(M) {
    this->P_M_2 = compute_power_fast(P, M - 2, M);
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
        (static_cast<long long>(elem) * this->computePreviousPowerP(this->nextPowerP)) + M) % M
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
    return this->P;
  }
  
  int getM() const {
    return this->M;
  }
  
  int getCurrentPowerP() const {
    return this->currentPowerP;
  }
  
private:
  int P;
  int M;
  int P_M_2; // P^(M-2) MOD M
  int nextPowerP;
  int currentPowerP;
  
  [[nodiscard]] int computeHash(int currentHash, const int &elem) const {
    return static_cast<int>(
      (static_cast<long long>(this->currentHash) * P + elem) % M
    );
  }
  
  [[nodiscard]] int computeNextPowerP(int currentPowerP) const {
    return static_cast<int>(
      (static_cast<long long>(currentPowerP) * P) % M
    );
  }
  
  [[nodiscard]] int computePreviousPowerP(int currentPowerP) const {
    return static_cast<int>(
      (static_cast<long long>(currentPowerP) * this->P_M_2) % M
    );
  }
};

#endif //SECOND_POLYNOMIAL_SINGLE_HASHER_H
