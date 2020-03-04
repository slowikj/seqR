// [[Rcpp::plugins("c++17")]]
#include "utils.h"
#include <queue>

int compute_power_fast(unsigned int base, unsigned int power, unsigned int modulo) {
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

