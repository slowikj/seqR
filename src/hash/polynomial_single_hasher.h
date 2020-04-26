#ifndef SECOND_POLYNOMIAL_SINGLE_HASHER_H
#define SECOND_POLYNOMIAL_SINGLE_HASHER_H

#include "single_hasher.h"
#include "../utils.h"

struct PolynomialSingleHasherConfig {
    int P;
    int M;

    PolynomialSingleHasherConfig(int P, int M) :
            P(P), M(M) {
    }

    PolynomialSingleHasherConfig() = delete;

    PolynomialSingleHasherConfig(const PolynomialSingleHasherConfig &) = default;

    PolynomialSingleHasherConfig &operator=(const PolynomialSingleHasherConfig &) = default;

};

class PolynomialSingleHasher : public SingleHasher {
public:
    PolynomialSingleHasher(PolynomialSingleHasherConfig &&config) :
            config(config) {
        this->P_M_2 = computePowerFast(config.P, config.M - 2, config.M);
        this->initPowersP();
    }

    inline void append(const int &elem) override {
        this->currentHash = this->computeHash(this->currentHash, elem);
        this->currentPowerP = this->nextPowerP;
        this->nextPowerP = this->computeNextPowerP(this->nextPowerP);
    }

    inline void removeFirst(const int &elem) override {
        this->currentHash = static_cast<int>(
                (this->currentHash -
                 (static_cast<long long>(elem) * this->currentPowerP) + config.M) % config.M
        );
        this->nextPowerP = this->computePreviousPowerP(this->nextPowerP);
        this->currentPowerP = this->computePreviousPowerP(this->currentPowerP);
    }

    inline int getHash() const override {
        return SingleHasher::getHash();
    }

    inline int getHash(int position) const override {
        return this->computeHash(currentHash, position);
    }

    inline void clear() override {
        SingleHasher::clear();
        this->initPowersP();
    }

    inline int getCurrentPowerP() const {
        return this->currentPowerP;
    }

private:
    PolynomialSingleHasherConfig config;
    int P_M_2; // P^(M-2) MOD M
    int nextPowerP;
    int currentPowerP;

    inline int computeHash(int currentHash, const int &elem) const {
        return static_cast<int>(
                (static_cast<long long>(this->currentHash) * config.P + elem) % config.M
        );
    }

    inline int computeNextPowerP(int currentPowerP) const {
        return static_cast<int>(
                (static_cast<long long>(currentPowerP) * config.P) % config.M
        );
    }

    inline int computePreviousPowerP(int currentPowerP) const {
        if (currentPowerP == 1) {
            return 0;
        }
        return static_cast<int>(
                (static_cast<long long>(currentPowerP) * this->P_M_2) % config.M
        );
    }

    inline void initPowersP() {
        this->nextPowerP = 1;
        this->currentPowerP = 0;
    }
};

#endif
