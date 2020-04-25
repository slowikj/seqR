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
#include "app_conf.h"

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
        class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix count_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> &alphabetEncoding,
                                int sequencesNum,
                                SequenceGetter_t<input_vector_t> sequenceGetter,
                                int k,
                                bool positionalKMers,
                                InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter) {
    Rcpp::IntegerVector gaps(k - 1);
    auto parallelKMerCountingProc = [&k, &positionalKMers, &sequencesNum](
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> &encoding,
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

template<class alphabet_t,
        class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix count_kmers(alphabet_t &alphabet,
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
Rcpp::IntegerMatrix count_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                                       Rcpp::StringVector &alphabet,
                                       int k,
                                       bool positionalKMers) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    std::vector<std::string> alphabetVector;
    std::transform(std::begin(alphabet), std::end(alphabet), std::back_inserter(alphabetVector),
            [](auto& elem) -> std::string { return Rcpp::as<std::string>(elem); });
    return count_kmers<
            std::vector<std::string>,
            SafeMatrixSequenceWrapper<std::string>::Row,
            std::string,
            ENCODED_ELEM_T,
            std::hash<std::string>>(alphabetVector,
                                 sequenceMatrix.nrow(),
                                 getRcppMatrixRowGetter<std::string>(safeMatrixWrapper),
                                 k,
                                 positionalKMers,
                                 getStringToStringConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                        Rcpp::IntegerVector &alphabet,
                                        int k,
                                        bool positionalKMers) {
    SafeMatrixSequenceWrapper<int> safeMatrixWrapper(sequenceMatrix);
    std::vector<int> alphabetVector;
    std::transform(std::begin(alphabet), std::end(alphabet), std::back_inserter(alphabetVector),
            [](int elem) -> int { return elem; });
    return count_kmers<
            std::vector<int>,
            SafeMatrixSequenceWrapper<int>::Row,
            int,
            ENCODED_ELEM_T,
            std::hash<int>>(alphabetVector,
                            sequenceMatrix.nrow(),
                            getRcppMatrixRowGetter<int>(safeMatrixWrapper),
                            k,
                            positionalKMers,
                            getIntToStringConverter());
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                        Rcpp::NumericVector &alphabet,
                                        int k,
                                        bool positionalKMers) {
    SafeMatrixSequenceWrapper<double> safeMatrixWrapper(sequenceMatrix);
    std::vector<double> alphabetVector;
    std::transform(std::begin(alphabet), std::end(alphabet), std::back_inserter(alphabetVector),
            [](double elem) -> double { return elem; });
    return count_kmers<
            std::vector<double>,
            SafeMatrixSequenceWrapper<double>::Row,
            double,
            ENCODED_ELEM_T,
            std::hash<double>>(alphabetVector,
                               sequenceMatrix.nrow(),
                               getRcppMatrixRowGetter<double>(safeMatrixWrapper),
                               k,
                               positionalKMers,
                               getDoubleToStringConverter(3));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers_tidysq(Rcpp::List &sq,
                                       Rcpp::StringVector &alphabet,
                                       int k,
                                       bool positionalKMers) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, std::hash<unsigned char>>(
            alphabet,
            elementsEncoding
    ));
    auto encodedSequences = getEncodedTidysqSequences(sq);
    SafeTidysqSequencesWrapper safeSequencesWrapper(encodedSequences);

    return count_kmers<
            SafeTidysqSequencesWrapper::Row,
            unsigned char,
            unsigned char,
            std::hash<unsigned char>>(alphabetEncoding,
                                      sq.size(),
                                      getTidysqRowGetter(safeSequencesWrapper),
                                      k,
                                      positionalKMers,
                                      getEncodedTidySqItemToStringConverter(elementsEncoding));
}
