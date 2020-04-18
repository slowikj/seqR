#include <Rcpp.h>
#include <memory>
#include <functional>
#include <vector>
#include "kmer_counting_common.h"
#include "kmer_counter.h"
#include "input_to_string_item_converter.h"
#include "sequence_getter.h"
#include "tidysq_encoded_sequence.h"
#include "app_conf.h"

inline ComplexHasher createKMerComplexHasher() {
  std::vector<std::unique_ptr<SingleHasher>> singleHashers;
  singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(101, 1e9 + 33)));
  singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(97, 1e9 + 7)));
  ComplexHasher complexHasher(std::move(singleHashers));
  return complexHasher;
}

template <class input_vector_t,
          class input_elem_t,
          class encoded_elem_t,
          class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix count_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t>& alphabetEncoding,
                                int sequencesNum,
                                SequenceGetter_t<input_vector_t> sequenceGetter,
                                int k,
                                bool positionalKMers,
                                InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
  Rcpp::IntegerVector gaps(k-1);
  auto parallelKMerCountingProc = [&k, &positionalKMers, &sequencesNum](
    AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t>& encoding,
    SequenceGetter_t<input_vector_t> seqGetter
    ) -> std::vector<KMerCountsManager> {
      return std::move(
        parallelComputeKMerCounts<input_vector_t,
                                  input_elem_t,
                                  encoded_elem_t,
                                  alphabet_hasher_t>(
                                    k,
                                    positionalKMers,
                                    sequencesNum,
                                    seqGetter,
                                    encoding,
                                    []() -> ComplexHasher { return createKMerComplexHasher(); }));
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
Rcpp::IntegerMatrix count_kmers(alphabet_t& alphabet,
                                int sequencesNum,
                                SequenceGetter_t<input_vector_t> sequenceGetter,
                                int k,
                                bool positionalKMers,
                                InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
  auto alphabetEncoding = std::move(getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_hasher_t>(
    alphabet
  ));
  return std::move(count_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_hasher_t>(
      alphabetEncoding,
      sequencesNum,
      sequenceGetter,
      k,
      positionalKMers,
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
                Rcpp::StringVector::stored_type,
                ENCODED_ELEM_T,
                string_proxy_hasher>(
                                alphabet,
                                sequenceMatrix.nrow(),
                                getRcppMatrixRowGetter<Rcpp::StringMatrix, Rcpp::StringMatrix::Row>(sequenceMatrix),
                                k,
                                positionalKMers,
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
                ENCODED_ELEM_T,
                std::hash<int>>(alphabet,
                                sequenceMatrix.nrow(),
                                getRcppMatrixRowGetter<Rcpp::IntegerMatrix, Rcpp::IntegerMatrix::Row>(sequenceMatrix),
                                k,
                                positionalKMers,
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
                ENCODED_ELEM_T,
                std::hash<double>>(alphabet,
                                  sequenceMatrix.nrow(),
                                  getRcppMatrixRowGetter<Rcpp::NumericMatrix, Rcpp::NumericMatrix::Row>(sequenceMatrix),
                                  k,
                                  positionalKMers,
                                  getDoubleToStringConverter(3))
  );
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_tidysq(Rcpp::StringVector& alphabet,
                                       Rcpp::List& sq,
                                       int k,
                                       bool positionalKMers) {
  Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
  auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, std::hash<unsigned char>>(
    alphabet,
    elementsEncoding
  ));
  auto encodedSequences = getEncodedTidysqSequences(sq);
  return std::move(
    count_kmers<TidysqEncodedSequenceProxy,
                unsigned char,
                unsigned char,
                std::hash<unsigned char>>(alphabetEncoding,
                                          sq.size(),
                                          getTidysqRowGetter(encodedSequences),
                                          k,
                                          positionalKMers,
                                          getEncodedTidySqItemToStringConverter(elementsEncoding))
  );
}
