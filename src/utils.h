#ifndef SECOND_POWER_HELPER_H
#define SECOND_POWER_HELPER_H

#include<queue>

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

template<class T>
void clear(std::queue<T>& q) {
  std::queue<T> empty;
  std::swap(q, empty);
}

#endif //SECOND_POWER_HELPER_H
