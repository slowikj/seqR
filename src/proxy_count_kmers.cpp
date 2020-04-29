// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <functional>
#include <vector>
#include "kmer_counting_common.h"
#include "kmer_counter.h"
#include "input_to_string_item_converter.h"
#include "sequence_getter.h"
#include "tidysq_encoded_sequence.h"
#include "rcpp_to_cpp_converters.h"

inline ComplexHasher createKMerComplexHasher() {
    std::vector<std::unique_ptr<SingleHasher>> singleHashers;
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(101, 1e9 + 33)));
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(97, 1e9 + 7)));
    ComplexHasher complexHasher(std::move(singleHashers));
    return complexHasher;
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::IntegerMatrix count_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                int sequencesNum,
                                SequenceGetter_t<input_vector_t> sequenceGetter,
                                int k,
                                bool positionalKMers,
                                bool withKMerCounts,
                                InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
    std::vector<int> gaps(k - 1);
    auto parallelKMerCountingProc = [&k, &positionalKMers, &sequencesNum, withKMerCounts](
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding,
            SequenceGetter_t<input_vector_t> seqGetter
    ) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                parallelComputeKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                        k,
                        positionalKMers,
                        withKMerCounts,
                        sequencesNum,
                        seqGetter,
                        encoding,
                        []() -> ComplexHasher { return createKMerComplexHasher(); }));
    };
    return std::move(
            getKMerCountsMatrix<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
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
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::IntegerMatrix count_kmers(alphabet_t &alphabet,
                                int sequencesNum,
                                SequenceGetter_t<input_vector_t> sequenceGetter,
                                int k,
                                bool positionalKMers,
                                bool withKMerCounts,
                                InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabet
            ));
    return std::move(
            count_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    alphabetEncoding,
                    sequencesNum,
                    sequenceGetter,
                    k,
                    positionalKMers,
                    withKMerCounts,
                    inputToStringItemConverter
            ));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                                      Rcpp::StringVector &alphabet,
                                      int k,
                                      bool positionalKMers,
                                      bool withKMerCounts) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    std::vector<std::string> convertedAlphabet = std::move(
            convertRcppVector<std::string, Rcpp::StringVector>(alphabet));
    return count_kmers<
            std::vector<std::string>,
            SafeMatrixSequenceWrapper<std::string>::Row,
            std::string,
            short,
            UnorderedMapWrapper,
            UnorderedMapWrapper>(convertedAlphabet,
                                 sequenceMatrix.nrow(),
                                 getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
                                 k,
                                 positionalKMers,
                                 withKMerCounts,
                                 getStringToStringConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                       Rcpp::IntegerVector &alphabet,
                                       int k,
                                       bool positionalKMers,
                                       bool withKMerCounts) {
    std::vector<int> convertedAlphabet = std::move(convertRcppVector<int, Rcpp::IntegerVector>(alphabet));
    return count_kmers<
            std::vector<int>,
            RcppParallel::RMatrix<int>::Row,
            int,
            short,
            UnorderedMapWrapper,
            UnorderedMapWrapper>(convertedAlphabet,
                                 sequenceMatrix.nrow(),
                                 getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix),
                                 k,
                                 positionalKMers,
                                 withKMerCounts,
                                 getIntToStringConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                       Rcpp::NumericVector &alphabet,
                                       int k,
                                       bool positionalKMers,
                                       bool withKMerCounts) {
    std::vector<double> convertedAlphabet = std::move(convertRcppVector<double, Rcpp::NumericVector>(alphabet));
    return count_kmers<
            std::vector<double>,
            RcppParallel::RMatrix<double>::Row,
            double,
            short,
            UnorderedMapWrapper,
            UnorderedMapWrapper>(convertedAlphabet,
                                 sequenceMatrix.nrow(),
                                 getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix),
                                 k,
                                 positionalKMers,
                                 withKMerCounts,
                                 getDoubleToStringConverter(3));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_tidysq(Rcpp::List &sq,
                                      Rcpp::StringVector &alphabet,
                                      int k,
                                      bool positionalKMers,
                                      bool withKMerCounts) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, UnorderedMapWrapper>(
            alphabet,
            elementsEncoding
    ));
    auto encodedSequences = getEncodedTidysqSequences(sq);
    return count_kmers<
            RcppParallel::RVector<unsigned char>,
            unsigned char,
            unsigned char,
            UnorderedMapWrapper,
            UnorderedMapWrapper>(alphabetEncoding,
                                 sq.size(),
                                 getTidysqRVectorGetter(encodedSequences),
                                 k,
                                 positionalKMers,
                                 withKMerCounts,
                                 getEncodedTidySqItemToStringConverter(elementsEncoding));
}
