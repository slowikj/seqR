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

template <class input_matrix_t,
          class input_vector_t,
          class input_view_vector_t,
          class input_elem_t,
          class internal_elem_t,
          class encoded_elem_t>
Rcpp::IntegerMatrix count_gapped_kmers(input_vector_t& alphabet,
                                       input_matrix_t& sequenceMatrix,
                                       Rcpp::IntegerVector& gaps,
                                       bool positionalKMers) {
  auto alphabetEncoding = std::move(getEncoding(alphabet));
  std::function<std::vector<KMerCountsManager>()> parallelKMerCountingProc =
    [&gaps, &alphabetEncoding, &positionalKMers, &sequenceMatrix]
    () -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeGappedKMersCounts<input_matrix_t,
                                         input_view_vector_t,
                                         input_elem_t,
                                         internal_elem_t,
                                         encoded_elem_t>(
                                           gaps,
                                           positionalKMers,
                                           sequenceMatrix,
                                           alphabetEncoding));
    };
  return std::move(getKMerCountsMatrix<input_matrix_t>(
      sequenceMatrix,
      gaps,
      positionalKMers,
      parallelKMerCountingProc
  ));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_string(Rcpp::StringVector& alphabet,
                                              Rcpp::StringMatrix& sequenceMatrix,
                                              Rcpp::IntegerVector& gaps,
                                              bool positionalKMers) {
  return std::move(count_gapped_kmers<Rcpp::StringMatrix,
                                      Rcpp::StringVector,
                                      Rcpp::StringMatrix::Row,
                                      Rcpp::String::StringProxy,
                                      std::string,
                                      ENCODED_ELEM_T>(alphabet,
                                                      sequenceMatrix,
                                                      gaps,
                                                      positionalKMers));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_integer(Rcpp::IntegerVector& alphabet,
                                               Rcpp::IntegerMatrix& sequenceMatrix,
                                               Rcpp::IntegerVector& gaps,
                                               bool positionalKMers) {
  return std::move(count_gapped_kmers<Rcpp::IntegerMatrix,
                                      Rcpp::IntegerVector,
                                      Rcpp::IntegerMatrix::Row,
                                      int,
                                      int,
                                      ENCODED_ELEM_T>(alphabet,
                                                      sequenceMatrix,
                                                      gaps,
                                                      positionalKMers));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_numeric(Rcpp::NumericVector& alphabet,
                                               Rcpp::NumericMatrix& sequenceMatrix,
                                               Rcpp::IntegerVector& gaps,
                                               bool positionalKMers) {
  return std::move(count_gapped_kmers<Rcpp::NumericMatrix,
                                      Rcpp::NumericVector,
                                      Rcpp::NumericMatrix::Row,
                                      double,
                                      double,
                                      ENCODED_ELEM_T>(alphabet,
                                                      sequenceMatrix,
                                                      gaps,
                                                      positionalKMers));
}
