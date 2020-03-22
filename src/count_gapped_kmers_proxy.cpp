#include "gapped_kmer_counter.h"
#include "kmer_counts_manager.h"
#include <Rcpp.h>
#include <memory>

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

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers(Rcpp::StringVector& alphabet,
                                       Rcpp::StringMatrix& sequenceMatrix,
                                       Rcpp::IntegerVector& gaps,
                                       bool positionalKMers) {
  auto alphabetEncoding = std::move(getEncoding(alphabet));
  std::function<std::vector<KMerCountsManager>()> parallelKMerCountingProc =
    [&gaps, &alphabetEncoding, &positionalKMers, &sequenceMatrix]
    () -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeGappedKMersCounts<Rcpp::StringMatrix,
                                         Rcpp::StringMatrix::Row,
                                         Rcpp::String::StringProxy,
                                         std::string,
                                         ENCODED_ELEM_T>(
                                           gaps,
                                           positionalKMers,
                                           sequenceMatrix,
                                           alphabetEncoding));
    };
    return getKMerCountsMatrix<Rcpp::StringMatrix>(
      sequenceMatrix,
      gaps,
      positionalKMers,
      parallelKMerCountingProc
    );
}
