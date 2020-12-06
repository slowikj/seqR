#ifndef UTILS_H
#define UTILS_H

#include <queue>
#include <vector>
#include <algorithm>

namespace util {

    template<class T>
    inline void clear(std::queue<T> &q) {
        std::queue<T> empty;
        std::swap(q, empty);
    }

    inline unsigned long long
    computePowerFast(unsigned long long base, unsigned long long power, unsigned long long modulo) {
        long long res = 1;
        long long current_base_power = base;
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
}

#endif
