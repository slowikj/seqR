#include<Rcpp.h>
#include<memory>
#include "kmer_counter.h"
#include <functional>
#include <vector>
#include "kmer_counting_common.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers(Rcpp::StringVector& alphabet,
                                Rcpp::StringMatrix& sequenceMatrix,
                                int k,
                                bool positionalKMers) {
  Rcpp::IntegerVector gaps(k-1);
  auto alphabetEncoding = getEncoding(alphabet);
  std::function<std::vector<KMerCountsManager>()> parallelKMerCountingProc =
    [&k, &positionalKMers, &sequenceMatrix, &alphabetEncoding]
    () -> std::vector<KMerCountsManager> {
    return std::move(
      parallelComputeKMerCounts<Rcpp::StringMatrix,
                                Rcpp::StringMatrix::Row,
                                Rcpp::String::StringProxy,
                                std::string,
                                ENCODED_ELEM_T>(
                                  k,
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
