#pragma once

#include <algorithm>
#include <queue>
#include <vector>

namespace util {

template <class T>
inline void clear(std::queue<T> &q) {
  std::queue<T> empty;
  std::swap(q, empty);
}

template <class T>
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
    for (std::size_t v_i = 1; v_i < v.size(); ++v_i) {
      res[v_i] = res[v_i - 1] + v[v_i] + 1;
    }
  }
  return res;
}

inline int getSum(const std::vector<int> &v) {
  return std::accumulate(std::begin(v), std::end(v), 0);
}

inline std::size_t getKMerRange(const std::vector<int> &gaps) {
  return util::getSum(gaps) + gaps.size() + 1;  // Rcpp::sum(gaps + 1) + 1
}

inline std::size_t getIntervalLength(const std::pair<std::size_t, std::size_t> &interval) {
  return interval.second - interval.first + 1;
}

// encapsulated for further optimizations (e.g., a custom, optimized modulo operation)
class ModuloComputer {
 public:
  explicit ModuloComputer(uint64_t M) : M(M) {}

  inline uint64_t get(uint64_t a) const {
    return a % M;
  }

 private:
  uint64_t M;
};

}  // namespace util
