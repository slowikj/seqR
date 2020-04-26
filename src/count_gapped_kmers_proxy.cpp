// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <functional>
#include "gapped_kmer_counter.h"
#include "kmer_counting_common.h"
#include "input_to_string_item_converter.h"
#include "sequence_getter.h"
#include "app_conf.h"
#include "tidysq_encoded_sequence.h"
#include "rcpp_to_cpp_converters.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_contiguous_intervals_matrix(const Rcpp::IntegerVector &gaps) {
    auto intervals = getContiguousIntervals(gaps);
    Rcpp::IntegerMatrix res(intervals.size(), 2);
    for (int i = 0; i < intervals.size(); ++i) {
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

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix
count_gapped_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> &alphabetEncoding,
                   int sequencesNum,
                   SequenceGetter_t<input_vector_t> sequenceGetter,
                   const std::vector<int> &gaps,
                   bool positionalKMers,
                   InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
    auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
    auto parallelKMerCountingProc =
            [&gaps, &positionalKMers, &sequencesNum, &hasherConfigs](
                    AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> &encoding,
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

template<class alphabet_t,
        class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix count_gapped_kmers(alphabet_t &alphabet,
                                       int sequencesNum,
                                       SequenceGetter_t<input_vector_t> sequenceGetter,
                                       const std::vector<int> &gaps,
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
Rcpp::IntegerMatrix count_gapped_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                                              Rcpp::StringVector &alphabet,
                                              Rcpp::IntegerVector &gaps,
                                              bool positionalKMers) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    std::vector<std::string> convertedAlphabet = std::move(
            convertRcppVector<std::string, Rcpp::StringVector>(alphabet));
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    return count_gapped_kmers<
            std::vector<std::string>,
            SafeMatrixSequenceWrapper<std::string>::Row,
            std::string,
            ENCODED_ELEM_T,
            std::hash<std::string>>(convertedAlphabet,
                                    sequenceMatrix.nrow(),
                                    getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
                                    gapsConverted,
                                    positionalKMers,
                                    getStringToStringConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                               Rcpp::IntegerVector &alphabet,
                                               Rcpp::IntegerVector &gaps,
                                               bool positionalKMers) {
    std::vector<int> convertedAlphabet = std::move(convertRcppVector<int, Rcpp::IntegerVector>(alphabet));
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    return count_gapped_kmers<
            std::vector<int>,
            RcppParallel::RMatrix<int>::Row,
            int,
            ENCODED_ELEM_T,
            std::hash<int>>(convertedAlphabet,
                            sequenceMatrix.nrow(),
                            getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix),
                            gapsConverted,
                            positionalKMers,
                            getIntToStringConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                               Rcpp::NumericVector &alphabet,
                                               Rcpp::IntegerVector &gaps,
                                               bool positionalKMers) {
    std::vector<double> convertedAlphabet = std::move(convertRcppVector<double, Rcpp::NumericVector>(alphabet));
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    return count_gapped_kmers<
            std::vector<double>,
            RcppParallel::RMatrix<double>::Row,
            double,
            ENCODED_ELEM_T,
            std::hash<double>>(convertedAlphabet,
                               sequenceMatrix.nrow(),
                               getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix),
                               gapsConverted,
                               positionalKMers,
                               getDoubleToStringConverter(3));
}


//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_gapped_kmers_tidysq(Rcpp::List &sq,
                                              Rcpp::StringVector &alphabet,
                                              Rcpp::IntegerVector &gaps,
                                              bool positionalKMers) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, std::hash<unsigned char>>(
            alphabet,
            elementsEncoding
    ));
    auto encodedSequences = getEncodedTidysqSequences(sq);
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    return count_gapped_kmers<
            RcppParallel::RVector<unsigned char>,
            unsigned char,
            unsigned char,
            std::hash<unsigned char>>(alphabetEncoding,
                                      sq.size(),
                                      getTidysqRVectorGetter(encodedSequences),
                                      gapsConverted,
                                      positionalKMers,
                                      getEncodedTidySqItemToStringConverter(elementsEncoding));
}
