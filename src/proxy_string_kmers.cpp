#include <Rcpp.h>
#include <string>
#include <vector>
#include "kmer_counting_result.h"
#include "alphabet_encoder.h"
#include "dictionary/unordered_map_wrapper.h"
#include "rcpp_to_cpp_converters.h"
#include "safe_sequences_wrapper.h"
#include "kmer_task_config.h"
#include "config_common.h"
#include "config_kmer_counting.h"
#include "config_gapped_kmer_counting.h"
#include "kmer_task_solver.h"

class StringKMersComputer {
public:
    using ENCODED_ELEM_T = short;

    explicit StringKMersComputer(Rcpp::StringVector alphabet) {
        initAlphabetEncoding(alphabet);
    }

    void addKMers(Rcpp::StringMatrix sequenceMatrix,
                  int k,
                  bool positionalKMers,
                  bool withKMerCounts,
                  const std::string kmerDictionaryName) {
        std::vector<int> gaps(k - 1);
        std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
        addKMersCommon(sequenceMatrix, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, algorithmParams,
                       result);
    }

    void addGappedKMers(Rcpp::StringMatrix sequenceMatrix,
                        Rcpp::IntegerVector gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName) {
        auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
        auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
        addKMersCommon(sequenceMatrix, gapsConverted, positionalKMers, withKMerCounts, kmerDictionaryName,
                       hasherConfigs, result);
    }

    Rcpp::List toList() const {
        return result.toRcppList();
    }

private:
    KMerCountingResult result;
    AlphabetEncoding<std::string, ENCODED_ELEM_T, UnorderedMapWrapper> alphabetEncoding;

    inline void initAlphabetEncoding(Rcpp::StringVector &alphabet) {
        this->alphabetEncoding = std::move(
                prepareAlphabetEncodingFromRcpp<Rcpp::StringVector, std::string, ENCODED_ELEM_T, UnorderedMapWrapper>(
                        alphabet));
    }

    template<class algorithm_params_t>
    inline
    void addKMersCommon(Rcpp::StringMatrix &sequenceMatrix,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        KMerCountingResult &kMerCountingResult) {
        SafeMatrixSequenceWrapper<std::string> safeMatrixWrapper(sequenceMatrix);
        KMerTaskConfig<SafeMatrixSequenceWrapper<std::string>::Row, std::string> kMerTaskConfig(
                sequenceMatrix.nrow(),
                getSafeMatrixRowGetter<std::string>(safeMatrixWrapper),
                gaps,
                positionalKMers,
                withKMerCounts,
                getStringToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                SafeMatrixSequenceWrapper<std::string>::Row,
                std::string,
                ENCODED_ELEM_T,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    this->alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    }
};

RCPP_MODULE(string_kmers_computer) {
    Rcpp::class_<StringKMersComputer>("StringKMersComputer")
            .constructor<Rcpp::StringVector>()
            .method("countKMers", &StringKMersComputer::addKMers)
            .method("countGappedKMers", &StringKMersComputer::addGappedKMers)
            .method("toList", &StringKMersComputer::toList);
}
