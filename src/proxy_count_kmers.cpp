// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <vector>
#include "kmer_counting_common.h"
#include "kmer_counter.h"
#include "input_to_string_item_converter.h"
#include "sequence_getter.h"
#include "tidysq_encoded_sequence.h"
#include "rcpp_to_cpp_converters.h"
#include "dictionary/linear_list_dictionary.h"

inline ComplexHasher createKMerComplexHasher() {
    std::vector<std::unique_ptr<SingleHasher>> singleHashers;
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(101, 1e9 + 33)));
    singleHashers.emplace_back(new PolynomialSingleHasher(PolynomialSingleHasherConfig(97, 1e9 + 7)));
    ComplexHasher complexHasher(std::move(singleHashers));
    return complexHasher;
}

template<class vector_t, class elem_t>
struct KMerTaskConf {
    int sequencesNum;
    SequenceGetter_t<vector_t> sequenceGetter;
    int k;
    bool positionalKMers;
    bool withKMerCounts;
    InputToStringItemConverter_t<elem_t> inputToStringItemConverter;

    KMerTaskConf(int sequencesNum,
                 SequenceGetter_t<vector_t> sequenceGetter,
                 int k,
                 bool positionalKMers,
                 bool withKMerCounts,
                 InputToStringItemConverter_t<elem_t> inputToStringItemConverter) :
            sequencesNum(sequencesNum),
            sequenceGetter(sequenceGetter),
            k(k),
            positionalKMers(positionalKMers),
            withKMerCounts(withKMerCounts),
            inputToStringItemConverter(inputToStringItemConverter) {
    }
};

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t,
        template<typename key, typename value, typename...> class kmer_dictionary_t>
inline
Rcpp::IntegerMatrix count_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                KMerTaskConf<input_vector_t, input_elem_t> &kMerTaskConf) {
    std::vector<int> gaps(kMerTaskConf.k - 1);
    auto parallelKMerCountingProc = [&kMerTaskConf](
            AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding,
            SequenceGetter_t<input_vector_t> seqGetter
    ) -> std::vector<KMerManager<kmer_dictionary_t>> {
        return std::move(
                parallelComputeKMers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                        kMerTaskConf.k,
                        kMerTaskConf.positionalKMers,
                        kMerTaskConf.withKMerCounts,
                        kMerTaskConf.sequencesNum,
                        seqGetter,
                        encoding,
                        []() -> ComplexHasher { return createKMerComplexHasher(); }));
    };
    return std::move(
            getKMerCountsMatrix<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    alphabetEncoding,
                    kMerTaskConf.sequencesNum,
                    kMerTaskConf.sequenceGetter,
                    gaps,
                    kMerTaskConf.positionalKMers,
                    parallelKMerCountingProc,
                    kMerTaskConf.inputToStringItemConverter
            ));
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline
Rcpp::IntegerMatrix count_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                                KMerTaskConf<input_vector_t, input_elem_t> &kMerTaskConf,
                                const std::string &kmerDictionaryName) {
    if (kmerDictionaryName == "unordered_map") {
        return std::move(
                count_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, UnorderedMapWrapper>(
                        alphabetEncoding,
                        kMerTaskConf
                ));
    } else if (kmerDictionaryName == "linear_list") {
        return std::move(
                count_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, LinearListDictionary>(
                        alphabetEncoding,
                        kMerTaskConf
                ));
    } else {
        throw std::runtime_error("unsupported kmer dictionary name " + kmerDictionaryName);
    }
}

template<class alphabet_t,
        class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline
Rcpp::IntegerMatrix count_kmers(alphabet_t &alphabet,
                                KMerTaskConf<input_vector_t, input_elem_t> &kMerTaskConf,
                                const std::string &kmerDictionaryName) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabet
            ));
    return std::move(
            count_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabetEncoding,
                    kMerTaskConf,
                    kmerDictionaryName
            ));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                                      Rcpp::StringVector &alphabet,
                                      int k,
                                      bool positionalKMers,
                                      bool withKMerCounts,
                                      const std::string& kmerDictionaryName) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    std::vector<std::string> convertedAlphabet = std::move(
            convertRcppVector<std::string, Rcpp::StringVector>(alphabet));
    KMerTaskConf<SafeMatrixSequenceWrapper<std::string>::Row, std::string> kMerTaskConf(
            sequenceMatrix.nrow(),
            getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
            k,
            positionalKMers,
            withKMerCounts,
            getStringToStringConverter());
    return count_kmers<
            std::vector<std::string>,
            SafeMatrixSequenceWrapper<std::string>::Row,
            std::string,
            short,
            UnorderedMapWrapper>(convertedAlphabet,
                                 kMerTaskConf,
                                 kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                       Rcpp::IntegerVector &alphabet,
                                       int k,
                                       bool positionalKMers,
                                       bool withKMerCounts,
                                       const std::string& kmerDictionaryName) {
    std::vector<int> convertedAlphabet = std::move(convertRcppVector<int, Rcpp::IntegerVector>(alphabet));
    KMerTaskConf<RcppParallel::RMatrix<int>::Row, int> kMerTaskConf(
            sequenceMatrix.nrow(),
            getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix),
            k,
            positionalKMers,
            withKMerCounts,
            getIntToStringConverter());
    return count_kmers<
            std::vector<int>,
            RcppParallel::RMatrix<int>::Row,
            int,
            short,
            UnorderedMapWrapper>(convertedAlphabet,
                                 kMerTaskConf,
                                 kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                       Rcpp::NumericVector &alphabet,
                                       int k,
                                       bool positionalKMers,
                                       bool withKMerCounts,
                                       const std::string &kmerDictionaryName) {
    std::vector<double> convertedAlphabet = std::move(convertRcppVector<double, Rcpp::NumericVector>(alphabet));
    KMerTaskConf<RcppParallel::RMatrix<double>::Row, double> kMerTaskConf(
            sequenceMatrix.nrow(),
            getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix),
            k,
            positionalKMers,
            withKMerCounts,
            getDoubleToStringConverter(3));
    return count_kmers<
            std::vector<double>,
            RcppParallel::RMatrix<double>::Row,
            double,
            short,
            UnorderedMapWrapper>(convertedAlphabet,
                                 kMerTaskConf,
                                 kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_kmers_tidysq(Rcpp::List &sq,
                                      Rcpp::StringVector &alphabet,
                                      int k,
                                      bool positionalKMers,
                                      bool withKMerCounts,
                                      const std::string& kmerDictionaryName) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, UnorderedMapWrapper>(
            alphabet,
            elementsEncoding
    ));
    std::vector<std::string> safeElementsEncoding = convertRcppVector<std::string, Rcpp::StringVector>(
            elementsEncoding);
    auto encodedSequences = getEncodedTidysqSequences(sq);
    KMerTaskConf<RcppParallel::RVector<unsigned char>, unsigned char> kMerTaskConf(
            sq.size(),
            getTidysqRVectorGetter(encodedSequences),
            k,
            positionalKMers,
            withKMerCounts,
            getEncodedTidySqItemToStringConverter(safeElementsEncoding));
    return count_kmers<
            RcppParallel::RVector<unsigned char>,
            unsigned char,
            unsigned char,
            UnorderedMapWrapper>(alphabetEncoding,
                                 kMerTaskConf,
                                 kmerDictionaryName);
}
