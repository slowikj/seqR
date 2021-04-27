#pragma once

#include "single_hasher.h"
#include "../utils.h"
#include "globals.h"

namespace hashing
{

    struct PolynomialSingleHasherConfig
    {
        using elem_t = SingleHasher::hash_t;

        elem_t P, M;

        PolynomialSingleHasherConfig(elem_t P, elem_t M) : P(P), M(M)
        {
        }

        PolynomialSingleHasherConfig() = delete;

        PolynomialSingleHasherConfig(const PolynomialSingleHasherConfig &) = default;

        PolynomialSingleHasherConfig &operator=(const PolynomialSingleHasherConfig &) = default;
    };

    class PolynomialSingleHasher : public SingleHasher
    {
    public:
        explicit PolynomialSingleHasher(PolynomialSingleHasherConfig &&config) : config(config), moduloMComputer(config.M)
        {
            this->P_M_2 = util::computePowerFast(config.P, config.M - 2, config.M);
            this->initPowersP();
            this->M = config.M;
        }

        inline void append(elem_t elem) override
        {
            this->currentHash = this->computeHash(this->currentHash, elem);
            this->currentPowerP = this->nextPowerP;
            this->nextPowerP = this->computeNextPowerP(this->nextPowerP);
        }

        inline void removeFirst(elem_t elem) override
        {
            this->currentHash = this->currentHash + M - moduloMComputer.get(this->currentPowerP * elem);
            if (this->currentHash > config.M)
            {
                this->currentHash -= config.M;
            }
            this->nextPowerP = this->computePreviousPowerP(this->nextPowerP);
            this->currentPowerP = this->computePreviousPowerP(this->currentPowerP);
        }

        [[nodiscard]] inline hash_t getHash() const override
        {
            return SingleHasher::getHash();
        }

        inline void clear() override
        {
            SingleHasher::clear();
            this->initPowersP();
        }

        [[nodiscard]] inline hash_t getCurrentPowerP() const
        {
            return this->currentPowerP;
        }

    private:
        PolynomialSingleHasherConfig config;
        hash_t P_M_2; // P^(M-2) MOD M
        hash_t nextPowerP;
        hash_t currentPowerP;
        uint64_t M;
        util::ModuloComputer moduloMComputer;

        [[nodiscard]] inline hash_t computeHash(hash_t currentHash, elem_t elem) const
        {
            return moduloMComputer.get(this->currentHash * config.P + elem);
        }

        [[nodiscard]] inline hash_t computeNextPowerP(hash_t currentPowerP) const
        {
            return moduloMComputer.get(currentPowerP * config.P);
        }

        [[nodiscard]] inline hash_t computePreviousPowerP(hash_t currentPowerP) const
        {
            if (currentPowerP == 1)
            {
                return 0;
            }
            return moduloMComputer.get(currentPowerP * this->P_M_2);
        }

        inline void initPowersP()
        {
            this->nextPowerP = 1;
            this->currentPowerP = 0;
        }
    };
}
