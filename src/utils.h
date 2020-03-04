#ifndef SECOND_POWER_HELPER_H
#define SECOND_POWER_HELPER_H

#include<queue>

int compute_power_fast(unsigned int base, unsigned int power, unsigned int modulo);

template<class T>
void clear(std::queue<T>& q) {
  std::queue<T> empty;
  std::swap(q, empty);
}

#endif //SECOND_POWER_HELPER_H
