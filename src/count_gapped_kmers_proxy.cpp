#include "gapped_kmer_counter.h"
#include <Rcpp.h>

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_contiguous_intervals_matrix(const Rcpp::IntegerVector& gaps) {
  auto intervals = getContiguousIntervals(gaps);
  Rcpp::IntegerMatrix res(intervals.size(), 2);
  for(int i = 0; i < intervals.size(); ++i) {
    res(i, 0) = intervals[i].first + 1;
    res(i, 1) = intervals[i].second + 1;
  }
  return res;
}
