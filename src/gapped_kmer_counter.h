#ifndef GAPPED_KMER_COUNTER_H
#define GAPPED_KMER_COUNTER_H

#include "Rcpp.h"
#include <vector>
#include <utility>

std::vector<std::pair<int, int>> getContiguousIntervals(const Rcpp::IntegerVector& gaps);

#endif