#ifndef SECOND_POLYNOMIAL_SINGLE_HASHER_H
#define SECOND_POLYNOMIAL_SINGLE_HASHER_H

#include "single_hasher.h"
#include "../utils.h"
#include "globals.h"
#include "../../inst/include/fast_modulo.h"

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
            this->fastMod_M = config.M;
            this->fastMod_M_conv = fastmod::computeM_u64(this->fastMod_M);
        }

        inline void append(const int &elem) override {
            this->currentHash = this->computeHash(this->currentHash, elem);
            this->currentPowerP = this->nextPowerP;
            this->nextPowerP = this->computeNextPowerP(this->nextPowerP);
        }

        inline void removeFirst(const int &elem) override {
            this->currentHash = this->currentHash + fastMod_M - getModuloM(this->currentPowerP * elem);
            if (this->currentHash > config.M) {
                this->currentHash -= config.M;
            }
            this->nextPowerP = this->computePreviousPowerP(this->nextPowerP);
            this->currentPowerP = this->computePreviousPowerP(this->currentPowerP);
        }

        inline config::single_hash_t getHash() const override {
            return SingleHasher::getHash();
        }

        inline void clear() override {
            SingleHasher::clear();
            this->initPowersP();
        }

        inline config::single_hash_t getCurrentPowerP() const {
            return this->currentPowerP;
        }

    private:
        PolynomialSingleHasherConfig config;
        config::single_hash_t P_M_2; // P^(M-2) MOD M
        config::single_hash_t nextPowerP;
        config::single_hash_t currentPowerP;
        uint64_t fastMod_M;
        __uint128_t fastMod_M_conv;

        [[nodiscard]] inline config::single_hash_t
        computeHash(config::single_hash_t currentHash, const int &elem) const {
            return getModuloM(this->currentHash * config.P + elem);
        }

        [[nodiscard]] inline config::single_hash_t computeNextPowerP(config::single_hash_t currentPowerP) const {
            return getModuloM(currentPowerP * config.P);
        }

        [[nodiscard]] inline config::single_hash_t
        computePreviousPowerP(config::single_hash_t currentPowerP) const {
            if (currentPowerP == 1) {
                return 0;
            }
            return getModuloM(currentPowerP * this->P_M_2);
        }

        inline void initPowersP() {
            this->nextPowerP = 1;
            this->currentPowerP = 0;
        }

        inline uint64_t getModuloM(uint64_t a) const {
            return fastmod::fastmod_u64(a, fastMod_M_conv, fastMod_M);
        }
    };
}


#endif
