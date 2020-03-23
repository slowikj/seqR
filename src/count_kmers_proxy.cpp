#include <Rcpp.h>
#include <memory>
#include "kmer_counter.h"
#include <functional>
#include <vector>
#include "kmer_counting_common.h"

template <class input_matrix_t,
          class input_vector_t,
          class input_view_vector_t,
          class input_elem_t,
          class internal_elem_t,
          class encoded_elem_t>
Rcpp::IntegerMatrix count_kmers(input_vector_t& alphabet,
                                input_matrix_t& sequenceMatrix,
                                int k,
                                bool positionalKMers) {
  Rcpp::IntegerVector gaps(k-1);
  auto alphabetEncoding = getEncoding(alphabet);
  std::function<std::vector<KMerCountsManager>()> parallelKMerCountingProc =
    [&k, &positionalKMers, &sequenceMatrix, &alphabetEncoding]
    () -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeKMerCounts<input_matrix_t,
                                  input_view_vector_t,
                                  input_elem_t,
                                  internal_elem_t,
                                  ENCODED_ELEM_T>(
                                    k,
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
Rcpp::IntegerMatrix count_kmers_string(Rcpp::StringVector& alphabet,
                                       Rcpp::StringMatrix& sequenceMatrix,
                                       int k,
                                       bool positionalKMers) {
  return std::move(count_kmers<Rcpp::StringMatrix,
                               Rcpp::StringVector,
                               Rcpp::StringMatrix::Row,
                               Rcpp::String::StringProxy,
                               std::string,
                               ENCODED_ELEM_T>(alphabet,
                                               sequenceMatrix,
                                               k,
                                               positionalKMers)
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_integer(Rcpp::IntegerVector& alphabet,
                                        Rcpp::IntegerMatrix& sequenceMatrix,
                                        int k,
                                        bool positionalKMers) {
  return std::move(count_kmers<Rcpp::IntegerMatrix,
                               Rcpp::IntegerVector,
                               Rcpp::IntegerMatrix::Row,
                               int,
                               int,
                               ENCODED_ELEM_T>(alphabet,
                                               sequenceMatrix,
                                               k,
                                               positionalKMers)
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_numeric(Rcpp::NumericVector& alphabet,
                                        Rcpp::NumericMatrix& sequenceMatrix,
                                        int k,
                                        bool positionalKMers) {
  return std::move(count_kmers<Rcpp::NumericMatrix,
                               Rcpp::NumericVector,
                               Rcpp::NumericMatrix::Row,
                               double,
                               double,
                               ENCODED_ELEM_T>(alphabet,
                                               sequenceMatrix,
                                               k,
                                               positionalKMers)
  );
}
