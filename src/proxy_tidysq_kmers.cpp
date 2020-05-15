#include <Rcpp.h>
#include <string>
#include <vector>
#include "kmer_counting_result.h"
#include "alphabet_encoder.h"
#include "dictionary/unordered_map_wrapper.h"
#include "rcpp_to_cpp_converters.h"
#include "kmer_task_config.h"
#include "config_common.h"
#include "config_kmer_counting.h"
#include "config_gapped_kmer_counting.h"
#include "kmer_task_solver.h"
#include "tidysq_encoded_sequence.h"

class TidySqKMersComputer {
public:
    explicit TidySqKMersComputer(Rcpp::StringVector tidysqElementsEncoding,
                                 Rcpp::StringVector alphabet) {
        initAlphabetEncoding(tidysqElementsEncoding, alphabet);
        this->safeElementsEncoding = convertRcppVector<std::string, Rcpp::StringVector>(alphabet);
    }

    void addKMers(Rcpp::List sq,
                  int k,
                  bool positionalKMers,
                  bool withKMerCounts,
                  const std::string kmerDictionaryName) {
        std::vector<int> gaps(k - 1);
        std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
        addKMersCommon(sq, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, algorithmParams,
                       result);
    }

    void addGappedKMers(Rcpp::List sq,
                        Rcpp::IntegerVector &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string kmerDictionaryName) {
        auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
        auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
        addKMersCommon(sq, gapsConverted, positionalKMers, withKMerCounts, kmerDictionaryName,
                       hasherConfigs, result);
    }

    Rcpp::List toList() const {
        return result.toRcppList();
    }

private:
    AlphabetEncoding<unsigned char, unsigned char, UnorderedMapWrapper> alphabetEncoding;
    std::vector<std::string> safeElementsEncoding;
    KMerCountingResult result;

    inline
    void initAlphabetEncoding(Rcpp::StringVector &tidySqElementsEncoding,
                              Rcpp::StringVector &alphabet) {
        this->alphabetEncoding = std::move(prepareAlphabetEncodingForTidysq<unsigned char, UnorderedMapWrapper>(
                alphabet,
                tidySqElementsEncoding
        ));
    }

    template<class algorithm_params_t>
    inline
    void addKMersCommon(Rcpp::List &sq,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        KMerCountingResult &kMerCountingResult) {
        auto encodedSequences = getEncodedTidysqSequences(sq);
        KMerTaskConfig<RcppParallel::RVector<unsigned char>, unsigned char> kMerTaskConfig(
                sq.size(),
                getTidysqRVectorGetter(encodedSequences),
                gaps,
                positionalKMers,
                withKMerCounts,
                getEncodedTidySqItemToStringConverter(safeElementsEncoding),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                RcppParallel::RVector<unsigned char>,
                unsigned char,
                unsigned char,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    this->alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    }
};

RCPP_MODULE(tidySq_kmers_computer) {
    Rcpp::class_<TidySqKMersComputer>("TidySqKMersComputer")
            .constructor<Rcpp::StringVector, Rcpp::StringVector>()
            .method("addKMers", &TidySqKMersComputer::addKMers)
            .method("addGappedKMers", &TidySqKMersComputer::addGappedKMers)
            .method("toList", &TidySqKMersComputer::toList);
}
