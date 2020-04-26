#ifndef UTILS_H
#define UTILS_H

#include <queue>
#include <vector>
#include <algorithm>

template<class T>
inline void clear(std::queue<T> &q) {
    std::queue<T> empty;
    std::swap(q, empty);
}

inline int computePowerFast(unsigned int base, unsigned int power, unsigned int modulo) {
    long long res = 1;
    long long current_base_power = base;
    while (power > 0) {
        if (power & 1) {
            res = (res * current_base_power) % modulo;
        }
        power >>= 1;
        current_base_power = (current_base_power * current_base_power) % modulo;
    }
    return static_cast<int>(res);
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

#endif
