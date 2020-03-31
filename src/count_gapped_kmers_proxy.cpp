#include <Rcpp.h>
#include <memory>
#include "gapped_kmer_counter.h"
#include "kmer_counting_common.h"
#include "input_to_string_item_converter.h"
#include "input_to_internal_item_converter.h"
#include "sequence_getter.h"
#include "app_conf.h"

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

inline std::vector<PolynomialSingleHasherConfig> getGappedKMerHasherConfigs() {
  std::vector<PolynomialSingleHasherConfig> res;
  res.emplace_back(101, 1e9 + 7);
  res.emplace_back(97, 1e9 + 33);
  return std::move(res);
}

template <class alphabet_t,
          class input_vector_t,
          class input_elem_t,
          class internal_elem_t,
          class encoded_elem_t>
Rcpp::IntegerMatrix count_gapped_kmers(alphabet_t& alphabet,
                                       int sequencesNum,
                                       SequenceGetter_t<input_vector_t> sequenceGetter,
                                       Rcpp::IntegerVector& gaps,
                                       bool positionalKMers,
                                       InputToInternalItemConverter_t<input_elem_t, internal_elem_t> inputToInternalItemConverter,
                                       InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
  auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
  auto parallelKMerCountingProc =
    [&gaps, &positionalKMers, &sequencesNum, &hasherConfigs](
        AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& encoding,
        SequenceGetter_t<input_vector_t> seqGetter
    ) -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeGappedKMersCounts<input_vector_t,
                                         input_elem_t,
                                         internal_elem_t,
                                         encoded_elem_t>(
                                           gaps,
                                           positionalKMers,
                                           sequencesNum,
                                           seqGetter,
                                           encoding,
                                           hasherConfigs));
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
Rcpp::IntegerMatrix count_gapped_kmers_string(Rcpp::StringVector& alphabet,
                                              Rcpp::StringMatrix& sequenceMatrix,
                                              Rcpp::IntegerVector& gaps,
                                              bool positionalKMers) {
  return std::move(
    count_gapped_kmers<Rcpp::StringVector,
                       Rcpp::StringMatrix::Row,
                       Rcpp::String::StringProxy,
                       std::string,
                       ENCODED_ELEM_T>(alphabet,
                                       sequenceMatrix.nrow(),
                                       getRcppMatrixRowGetter<Rcpp::StringMatrix, Rcpp::StringMatrix::Row>(sequenceMatrix),
                                       gaps,
                                       positionalKMers,
                                       getRcppStringToStringConverter(),
                                       getRcppStringProxyToStringConverter())
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_integer(Rcpp::IntegerVector& alphabet,
                                               Rcpp::IntegerMatrix& sequenceMatrix,
                                               Rcpp::IntegerVector& gaps,
                                               bool positionalKMers) {
  return std::move(
    count_gapped_kmers<Rcpp::IntegerVector,
                       Rcpp::IntegerMatrix::Row,
                       int,
                       int,
                       ENCODED_ELEM_T>(alphabet,
                                       sequenceMatrix.nrow(),
                                       getRcppMatrixRowGetter<Rcpp::IntegerMatrix, Rcpp::IntegerMatrix::Row>(sequenceMatrix),
                                       gaps,
                                       positionalKMers,
                                       getIntToIntConverter(),
                                       getIntToStringConverter())
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_numeric(Rcpp::NumericVector& alphabet,
                                               Rcpp::NumericMatrix& sequenceMatrix,
                                               Rcpp::IntegerVector& gaps,
                                               bool positionalKMers) {
  return std::move(
    count_gapped_kmers<Rcpp::NumericVector,
                       Rcpp::NumericMatrix::Row,
                       double,
                       double,
                       ENCODED_ELEM_T>(alphabet,
                                       sequenceMatrix.nrow(),
                                       getRcppMatrixRowGetter<Rcpp::NumericMatrix, Rcpp::NumericMatrix::Row>(sequenceMatrix),
                                       gaps,
                                       positionalKMers,
                                       getDoubleToDoubleConverter(),
                                       getDoubleToStringConverter(3))
  );
}
