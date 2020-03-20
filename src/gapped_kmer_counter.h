#ifndef GAPPED_KMER_COUNTER_H
#define GAPPED_KMER_COUNTER_H

#include "Rcpp.h"
#include "alphabet_encoder.h"
#include <vector>
#include <memory>
#include <utility>

std::vector<std::pair<int, int>> getContiguousIntervals(const Rcpp::IntegerVector& gaps);

template<class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<int> prepareNotAllowedItemsPrefixCount(
    input_vector_t& sequence,
    const AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding) {
  std::vector<int> res;
  res.reserve(sequence.size());
  for(int i = 0; i < sequence.size(); ++i) {
    bool isNotPresent = !alphabetEncoding.isPresent(sequence[i]);
    res[i] = (i == 0) ?
      isNotPresent :
      res[i - 1] + isNotPresent;
  }
  return std::move(res);
}

#endif