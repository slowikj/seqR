#include <Rcpp.h>
#include <memory>
#include <functional>
#include "gapped_kmer_counter.h"
#include "kmer_counting_common.h"
#include "input_to_string_item_converter.h"
#include "hash/custom_hashers.h"
#include "sequence_getter.h"
#include "app_conf.h"
#include "tidysq_encoded_sequence.h"

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

template <class input_vector_t,
          class input_elem_t,
          class encoded_elem_t,
          class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix count_gapped_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t>& alphabetEncoding,
                                       int sequencesNum,
                                       SequenceGetter_t<input_vector_t> sequenceGetter,
                                       Rcpp::IntegerVector& gaps,
                                       bool positionalKMers,
                                       InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
  auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
  auto parallelKMerCountingProc =
    [&gaps, &positionalKMers, &sequencesNum, &hasherConfigs](
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t>& encoding,
        SequenceGetter_t<input_vector_t> seqGetter
    ) -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeGappedKMersCounts<input_vector_t, input_elem_t, encoded_elem_t, alphabet_hasher_t>(
          gaps,
          positionalKMers,
          sequencesNum,
         seqGetter,
          encoding,
          hasherConfigs
      ));
  };
  return std::move(getKMerCountsMatrix<input_vector_t, input_elem_t, encoded_elem_t, alphabet_hasher_t>(
      alphabetEncoding,
      sequencesNum,
      sequenceGetter,
      gaps,
      positionalKMers,
      parallelKMerCountingProc,
      inputToStringItemConverter
  ));
}

template <class alphabet_t,
          class input_vector_t,
          class input_elem_t,
          class encoded_elem_t,
          class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix count_gapped_kmers(alphabet_t& alphabet,
                                       int sequencesNum,
                                       SequenceGetter_t<input_vector_t> sequenceGetter,
                                       Rcpp::IntegerVector& gaps,
                                       bool positionalKMers,
                                       InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
  auto alphabetEncoding = std::move(getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_hasher_t>(
    alphabet
  ));
  return std::move(count_gapped_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_hasher_t>(
    alphabetEncoding,
    sequencesNum,
    sequenceGetter,
    gaps,
    positionalKMers,
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
                       Rcpp::StringVector::stored_type,
                       ENCODED_ELEM_T,
                       string_proxy_hasher>(alphabet,
                                            sequenceMatrix.nrow(),
                                            getRcppMatrixRowGetter<Rcpp::StringMatrix, Rcpp::StringMatrix::Row>(sequenceMatrix),
                                            gaps,
                                            positionalKMers,
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
                       ENCODED_ELEM_T,
                       std::hash<int>>(alphabet,
                                       sequenceMatrix.nrow(),
                                       getRcppMatrixRowGetter<Rcpp::IntegerMatrix, Rcpp::IntegerMatrix::Row>(sequenceMatrix),
                                       gaps,
                                       positionalKMers,
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
                       ENCODED_ELEM_T,
                       std::hash<double>>(alphabet,
                                          sequenceMatrix.nrow(),
                                          getRcppMatrixRowGetter<Rcpp::NumericMatrix, Rcpp::NumericMatrix::Row>(sequenceMatrix),
                                          gaps,
                                          positionalKMers,
                                          getDoubleToStringConverter(3))
  );
}


//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_tidysq(Rcpp::StringVector& alphabet,
                                              Rcpp::List& sq,
                                              Rcpp::IntegerVector& gaps,
                                              bool positionalKMers) {
  Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
  auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, std::hash<unsigned char>>(
    alphabet,
    elementsEncoding 
  ));
  auto encodedSequences = getEncodedTidysqSequences(sq);
  
  return std::move(
    count_gapped_kmers<Rcpp::RawVector,
                       unsigned char,
                       unsigned char,
                       std::hash<unsigned char>>(alphabetEncoding,
                                                 sq.size(),
                                                 getTidysqRowGetter(encodedSequences),
                                                 gaps,
                                                 positionalKMers,
                                                 getEncodedTidySqItemToStringConverter(elementsEncoding))
  );
}
