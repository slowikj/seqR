#include <Rcpp.h>
#include <memory>
#include "kmer_counter.h"
#include <functional>
#include <vector>
#include "kmer_counting_common.h"
#include "input_to_string_item_converter.h"
#include "input_to_internal_item_converter.h"
#include "sequence_getter.h"
#include "app_conf.h"

template <class alphabet_t,
          class input_vector_t,
          class input_elem_t,
          class internal_elem_t,
          class encoded_elem_t>
Rcpp::IntegerMatrix count_kmers(alphabet_t& alphabet,
                                int sequencesNum,
                                SequenceGetter_t<input_vector_t> sequenceGetter,
                                int k,
                                bool positionalKMers,
                                InputToInternalItemConverter_t<input_elem_t, internal_elem_t> inputToInternalItemConverter,
                                InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
  Rcpp::IntegerVector gaps(k-1);
  auto parallelKMerCountingProc = [&k, &positionalKMers, &sequencesNum](
    AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& encoding,
    SequenceGetter_t<input_vector_t> seqGetter
    ) -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeKMerCounts<input_vector_t,
                                  input_elem_t,
                                  internal_elem_t,
                                  ENCODED_ELEM_T>(
                                    k,
                                    positionalKMers,
                                    sequencesNum,
                                    seqGetter,
                                    encoding));
    };
  return std::move(getKMerCountsMatrix<alphabet_t, input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
      alphabet,
      sequencesNum,
      sequenceGetter,
      gaps,
      positionalKMers,
      parallelKMerCountingProc,
      inputToInternalItemConverter,
      inputToStringItemConverter
  ));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_string(Rcpp::StringVector& alphabet,
                                       Rcpp::StringMatrix& sequenceMatrix,
                                       int k,
                                       bool positionalKMers) {
  return std::move(
    count_kmers<Rcpp::StringVector,
                Rcpp::StringMatrix::Row,
                Rcpp::String::StringProxy,
                std::string,
                ENCODED_ELEM_T>(alphabet,
                                sequenceMatrix.nrow(),
                                getRcppMatrixRowGetter<Rcpp::StringMatrix, Rcpp::StringMatrix::Row>(sequenceMatrix),
                                k,
                                positionalKMers,
                                getRcppStringToStringConverter(),
                                getRcppStringProxyToStringConverter())
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_integer(Rcpp::IntegerVector& alphabet,
                                        Rcpp::IntegerMatrix& sequenceMatrix,
                                        int k,
                                        bool positionalKMers) {
  return std::move(
    count_kmers<Rcpp::IntegerVector,
                Rcpp::IntegerMatrix::Row,
                int,
                int,
                ENCODED_ELEM_T>(alphabet,
                                sequenceMatrix.nrow(),
                                getRcppMatrixRowGetter<Rcpp::IntegerMatrix, Rcpp::IntegerMatrix::Row>(sequenceMatrix),
                                k,
                                positionalKMers,
                                getIntToIntConverter(),
                                getIntToStringConverter())
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_numeric(Rcpp::NumericVector& alphabet,
                                        Rcpp::NumericMatrix& sequenceMatrix,
                                        int k,
                                        bool positionalKMers) {
  return std::move(
    count_kmers<Rcpp::NumericVector,
                Rcpp::NumericMatrix::Row,
                double,
                double,
                ENCODED_ELEM_T>(alphabet,
                                sequenceMatrix.nrow(),
                                getRcppMatrixRowGetter<Rcpp::NumericMatrix, Rcpp::NumericMatrix::Row>(sequenceMatrix),
                                k,
                                positionalKMers,
                                getDoubleToDoubleConverter(),
                                getDoubleToStringConverter(3))
  );
}
