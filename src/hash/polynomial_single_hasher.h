#ifndef SECOND_POLYNOMIAL_SINGLE_HASHER_H
#define SECOND_POLYNOMIAL_SINGLE_HASHER_H

#include "single_hasher.h"
#include "../utils.h"
#include "types.h"

namespace hashing {

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
            this->P_M_2 = util::computePowerFast(config.P, config.M - 2, config.M);
            this->initPowersP();
        }

        inline void append(const int &elem) override {
            this->currentHash = this->computeHash(this->currentHash, elem);
            this->currentPowerP = this->nextPowerP;
            this->nextPowerP = this->computeNextPowerP(this->nextPowerP);
        }

        inline void removeFirst(const int &elem) override {
            this->currentHash = ((this->currentHash -
                     (this->currentPowerP * elem) % config.M) + config.M) % config.M;
            this->nextPowerP = this->computePreviousPowerP(this->nextPowerP);
            this->currentPowerP = this->computePreviousPowerP(this->currentPowerP);
        }

        inline single_hash_t getHash() const override {
            return SingleHasher::getHash();
        }

        inline void clear() override {
            SingleHasher::clear();
            this->initPowersP();
        }

        inline single_hash_t getCurrentPowerP() const {
            return this->currentPowerP;
        }

    private:
        PolynomialSingleHasherConfig config;
        single_hash_t P_M_2; // P^(M-2) MOD M
        single_hash_t nextPowerP;
        single_hash_t currentPowerP;

        [[nodiscard]] inline single_hash_t computeHash(single_hash_t currentHash, const int &elem) const {
            return (this->currentHash * config.P + elem) % config.M;
        }

        [[nodiscard]] inline single_hash_t computeNextPowerP(single_hash_t currentPowerP) const {
            return (currentPowerP * config.P) % config.M;
        }

        [[nodiscard]] inline single_hash_t computePreviousPowerP(single_hash_t currentPowerP) const {
            if (currentPowerP == 1) {
                return 0;
            }
            return (currentPowerP * this->P_M_2) % config.M;
        }

        inline void initPowersP() {
            this->nextPowerP = 1;
            this->currentPowerP = 0;
        }
    };
}


#endif
