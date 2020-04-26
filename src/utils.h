#ifndef UTILS_H
#define UTILS_H

#include <queue>
#include <Rcpp.h>

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

inline Rcpp::IntegerVector getGapsAccumulated(const Rcpp::IntegerVector &gaps) {
    return static_cast<Rcpp::IntegerVector>(Rcpp::cumsum(gaps + 1));
}

#endif
