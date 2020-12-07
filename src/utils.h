#ifndef UTILS_H
#define UTILS_H

#include <queue>
#include <vector>
#include <algorithm>
#include "../inst/include/fast_modulo.h"

namespace util {

    template<class T>
    inline void clear(std::queue<T> &q) {
        std::queue<T> empty;
        std::swap(q, empty);
    }

    template<class T>
    inline T computePowerFast(T base, T power, T modulo) {
        T res = 1;
        T current_base_power = base;
        while (power > 0) {
            if (power & 1) {
                res = (res * current_base_power) % modulo;
            }
            power >>= 1;
            current_base_power = (current_base_power * current_base_power) % modulo;
        }
        return res;
    }

    inline std::vector<int> getGapsAccumulated(const std::vector<int> &v) {
        std::vector<int> res(v.size());
        if (v.size() > 0) {
            res[0] = v[0] + 1;
            for (int v_i = 1; v_i < v.size(); ++v_i) {
                res[v_i] = res[v_i - 1] + v[v_i] + 1;
            }
        }
        return std::move(res);
    }

    inline int getSum(const std::vector<int> &v) {
        return std::accumulate(std::begin(v), std::end(v), 0);
    }

    inline std::size_t getKMerRange(const std::vector<int> &gaps) {
        return util::getSum(gaps) + gaps.size() + 1; // Rcpp::sum(gaps + 1) + 1
    }

    inline int getIntervalLength(const std::pair<int, int> &interval) {
        return interval.second - interval.first + 1;
    }

#ifdef __SIZEOF_INT128__

    class ModuloComputer {
    public:
        explicit ModuloComputer(uint64_t M) :
                M(M),
                fastMod_M_conv(fastmod::computeM_u64(M)) {
        }

        inline uint64_t get(uint64_t a) const {
            return fastmod::fastmod_u64(a, fastMod_M_conv, M);
        }

    private:
        uint64_t M;
        __uint128_t fastMod_M_conv;
    };

#else
    class ModuloComputer {
    public:
        explicit MoModuloComputer(uintuint64_t M) : M(M) {}

        inline uint64_t get(uintuint64_t a) const {
            return a % M;
        }
    };
#endif


}

#endif
