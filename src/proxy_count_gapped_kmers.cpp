// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <functional>
#include "gapped_kmer_counter.h"
#include "kmer_counting_common.h"
#include "input_to_string_item_converter.h"
#include "sequence_getter.h"
#include "tidysq_encoded_sequence.h"
#include "rcpp_to_cpp_converters.h"
#include "dictionary/linear_list_dictionary.h"

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

template<class vector_t, class elem_t>
struct GappedKMerTaskConf {
    int sequencesNum;
    SequenceGetter_t<vector_t> sequenceGetter;
    const std::vector<int> &gaps;
    bool positionalKMers;
    bool withKMerCounts;
    InputToStringItemConverter_t<elem_t> inputToStringItemConverter;

    GappedKMerTaskConf(int sequencesNum,
                       SequenceGetter_t<vector_t> sequenceGetter,
                       const std::vector<int> &gaps,
                       bool positionalKMers,
                       bool withKMerCounts,
                       InputToStringItemConverter_t<elem_t> inputToStringItemConverter) :
            sequencesNum(sequencesNum),
            sequenceGetter(sequenceGetter),
            gaps(gaps),
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
Rcpp::IntegerMatrix
count_gapped_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   GappedKMerTaskConf<input_vector_t, input_elem_t> &taskConf) {
    auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
    auto parallelKMerCountingProc =
            [&taskConf, &hasherConfigs](
                    AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &encoding,
                    SequenceGetter_t<input_vector_t> seqGetter
            ) -> std::vector<KMerManager<kmer_dictionary_t>> {
                return std::move(
                        parallelComputeGappedKMersCounts<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                                taskConf.gaps,
                                taskConf.positionalKMers,
                                taskConf.withKMerCounts,
                                taskConf.sequencesNum,
                                seqGetter,
                                encoding,
                                hasherConfigs
                        ));
            };
    return std::move(
            getKMerCountsMatrix<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, kmer_dictionary_t>(
                    alphabetEncoding,
                    taskConf.sequencesNum,
                    taskConf.sequenceGetter,
                    taskConf.gaps,
                    taskConf.positionalKMers,
                    parallelKMerCountingProc,
                    taskConf.inputToStringItemConverter
            ));
}

template<class input_vector_t,
        class input_elem_t,
        class encoded_elem_t,
        template<typename input_t, typename encoded_t, typename...> class alphabet_dictionary_t>
inline
Rcpp::IntegerMatrix
count_gapped_kmers(AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_dictionary_t> &alphabetEncoding,
                   GappedKMerTaskConf<input_vector_t, input_elem_t> &taskConf,
                   const std::string &kmerDictionaryName) {
    if (kmerDictionaryName == "unordered_map") {
        return std::move(
                count_gapped_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, UnorderedMapWrapper>(
                        alphabetEncoding,
                        taskConf
                ));
    } else if (kmerDictionaryName == "linear_list") {
        return std::move(
                count_gapped_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t, LinearListDictionary>(
                        alphabetEncoding,
                        taskConf
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
Rcpp::IntegerMatrix count_gapped_kmers(alphabet_t &alphabet,
                                       GappedKMerTaskConf<input_vector_t, input_elem_t> &taskConf,
                                       const std::string &kmerDictionaryName) {
    auto alphabetEncoding = std::move(
            getAlphabetEncoding<alphabet_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabet
            ));
    return std::move(
            count_gapped_kmers<input_vector_t, input_elem_t, encoded_elem_t, alphabet_dictionary_t>(
                    alphabetEncoding,
                    taskConf,
                    kmerDictionaryName
            ));
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                                             Rcpp::StringVector &alphabet,
                                             Rcpp::IntegerVector &gaps,
                                             bool positionalKMers,
                                             bool withKMerCounts,
                                             const std::string &kmerDictionaryName) {
    SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
    std::vector<std::string> convertedAlphabet = std::move(
            convertRcppVector<std::string, Rcpp::StringVector>(alphabet));
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    GappedKMerTaskConf<SafeMatrixSequenceWrapper<std::string>::Row, std::string> taskConf(
            sequenceMatrix.nrow(),
            getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
            gapsConverted,
            positionalKMers,
            withKMerCounts,
            getStringToStringConverter()
    );
    return count_gapped_kmers<
            std::vector<std::string>,
            SafeMatrixSequenceWrapper<std::string>::Row,
            std::string,
            short,
            UnorderedMapWrapper>(convertedAlphabet,
                                 taskConf,
                                 kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                              Rcpp::IntegerVector &alphabet,
                                              Rcpp::IntegerVector &gaps,
                                              bool positionalKMers,
                                              bool withKMerCounts,
                                              const std::string &kmerDictionaryName) {
    std::vector<int> convertedAlphabet = std::move(convertRcppVector<int, Rcpp::IntegerVector>(alphabet));
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    GappedKMerTaskConf<RcppParallel::RMatrix<int>::Row, int> taskConf(
            sequenceMatrix.nrow(),
            getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix),
            gapsConverted,
            positionalKMers,
            withKMerCounts,
            getIntToStringConverter()
    );
    return count_gapped_kmers<
            std::vector<int>,
            RcppParallel::RMatrix<int>::Row,
            int,
            short,
            UnorderedMapWrapper>(convertedAlphabet,
                                 taskConf,
                                 kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                              Rcpp::NumericVector &alphabet,
                                              Rcpp::IntegerVector &gaps,
                                              bool positionalKMers,
                                              bool withKMerCounts,
                                              const std::string &kmerDictionaryName) {
    std::vector<double> convertedAlphabet = std::move(convertRcppVector<double, Rcpp::NumericVector>(alphabet));
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    GappedKMerTaskConf<RcppParallel::RMatrix<double>::Row, double> taskConf(
            sequenceMatrix.nrow(),
            getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix),
            gapsConverted,
            positionalKMers,
            withKMerCounts,
            getDoubleToStringConverter(3)
    );
    return count_gapped_kmers<
            std::vector<double>,
            RcppParallel::RMatrix<double>::Row,
            double,
            short,
            UnorderedMapWrapper>(convertedAlphabet,
                                 taskConf,
                                 kmerDictionaryName);
}


//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_tidysq(Rcpp::List &sq,
                                             Rcpp::StringVector &alphabet,
                                             Rcpp::IntegerVector &gaps,
                                             bool positionalKMers,
                                             bool withKMerCounts,
                                             const std::string &kmerDictionaryName) {
    Rcpp::StringVector elementsEncoding = sq.attr("alphabet");
    auto alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, UnorderedMapWrapper>(
            alphabet,
            elementsEncoding
    ));
    std::vector<std::string> safeElementsEncoding = convertRcppVector<std::string, Rcpp::StringVector>(
            elementsEncoding);
    auto encodedSequences = getEncodedTidysqSequences(sq);
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    GappedKMerTaskConf<RcppParallel::RVector<unsigned char>, unsigned char> taskConf(
            sq.size(),
            getTidysqRVectorGetter(encodedSequences),
            gapsConverted,
            positionalKMers,
            withKMerCounts,
            getEncodedTidySqItemToStringConverter(safeElementsEncoding)
    );
    return count_gapped_kmers<
            RcppParallel::RVector<unsigned char>,
            unsigned char,
            unsigned char,
            UnorderedMapWrapper>(alphabetEncoding,
                                 taskConf,
                                 kmerDictionaryName);
}
