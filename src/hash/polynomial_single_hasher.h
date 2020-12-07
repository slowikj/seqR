#ifndef SECOND_POLYNOMIAL_SINGLE_HASHER_H
#define SECOND_POLYNOMIAL_SINGLE_HASHER_H

#include "single_hasher.h"
#include "../utils.h"
#include "globals.h"
#include "../../inst/include/fast_modulo.h"

namespace hashing {

    struct PolynomialSingleHasherConfig {
        using elem_t = SingleHasher::hash_t;

        elem_t P, M;

        PolynomialSingleHasherConfig(elem_t P, elem_t M) :
                P(P), M(M) {
        }

        PolynomialSingleHasherConfig() = delete;

        PolynomialSingleHasherConfig(const PolynomialSingleHasherConfig &) = default;

        PolynomialSingleHasherConfig &operator=(const PolynomialSingleHasherConfig &) = default;
    };

    class PolynomialSingleHasher : public SingleHasher {
    public:
        explicit PolynomialSingleHasher(PolynomialSingleHasherConfig &&config):
                config(config) {
            this->P_M_2 = util::computePowerFast(config.P, config.M - 2, config.M);
            this->initPowersP();
            this->fastMod_M = config.M;
            this->fastMod_M_conv = fastmod::computeM_u64(this->fastMod_M);
        }

        inline void append(elem_t elem) override {
            this->currentHash = this->computeHash(this->currentHash, elem);
            this->currentPowerP = this->nextPowerP;
            this->nextPowerP = this->computeNextPowerP(this->nextPowerP);
        }

        inline void removeFirst(elem_t elem) override {
            this->currentHash = this->currentHash + fastMod_M - getModuloM(this->currentPowerP * elem);
            if (this->currentHash > config.M) {
                this->currentHash -= config.M;
            }
            this->nextPowerP = this->computePreviousPowerP(this->nextPowerP);
            this->currentPowerP = this->computePreviousPowerP(this->currentPowerP);
        }

        [[nodiscard]] inline hash_t getHash() const override {
            return SingleHasher::getHash();
        }

        inline void clear() override {
            SingleHasher::clear();
            this->initPowersP();
        }

        [[nodiscard]] inline hash_t getCurrentPowerP() const {
            return this->currentPowerP;
        }

    private:
        PolynomialSingleHasherConfig config;
        hash_t P_M_2; // P^(M-2) MOD M
        hash_t nextPowerP;
        hash_t currentPowerP;
        uint64_t fastMod_M;
        __uint128_t fastMod_M_conv;

        [[nodiscard]] inline hash_t computeHash(hash_t currentHash, elem_t elem) const {
            return getModuloM(this->currentHash * config.P + elem);
        }

        [[nodiscard]] inline hash_t computeNextPowerP(hash_t currentPowerP) const {
            return getModuloM(currentPowerP * config.P);
        }

        [[nodiscard]] inline hash_t computePreviousPowerP(hash_t currentPowerP) const {
            if (currentPowerP == 1) {
                return 0;
            }
            return getModuloM(currentPowerP * this->P_M_2);
        }

        inline void initPowersP() {
            this->nextPowerP = 1;
            this->currentPowerP = 0;
        }

        [[nodiscard]] inline uint64_t getModuloM(uint64_t a) const {
            return fastmod::fastmod_u64(a, fastMod_M_conv, fastMod_M);
        }
    };
}


#endif
