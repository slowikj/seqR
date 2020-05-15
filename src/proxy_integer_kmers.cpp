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

class IntegerKMersComputer {
public:
    using ENCODED_ELEM_T = short;

    explicit IntegerKMersComputer(Rcpp::IntegerVector alphabet) {
        initAlphabetEncoding(alphabet);
    }

    void addKMers(Rcpp::IntegerMatrix sequenceMatrix,
                  int k,
                  bool positionalKMers,
                  bool withKMerCounts,
                  const std::string kmerDictionaryName) {
        std::vector<int> gaps(k - 1);
        std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
        addKMersCommon(sequenceMatrix, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, algorithmParams,
                       result);
    }

    void addGappedKMers(Rcpp::IntegerMatrix sequenceMatrix,
                        Rcpp::IntegerVector gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string kmerDictionaryName) {
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
    AlphabetEncoding<int, ENCODED_ELEM_T, UnorderedMapWrapper> alphabetEncoding;

    inline void initAlphabetEncoding(Rcpp::IntegerVector &alphabet) {
        this->alphabetEncoding = std::move(
                prepareAlphabetEncodingFromRcpp<Rcpp::IntegerVector, int, ENCODED_ELEM_T, UnorderedMapWrapper>(
                        alphabet));
    }

    template<class algorithm_params_t>
    inline
    void addKMersCommon(Rcpp::IntegerMatrix &sequenceMatrix,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        KMerCountingResult &kMerCountingResult) {
        KMerTaskConfig<RcppParallel::RMatrix<int>::Row, int> kMerTaskConfig(
                sequenceMatrix.nrow(),
                getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix),
                gaps,
                positionalKMers,
                withKMerCounts,
                getIntToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                RcppParallel::RMatrix<int>::Row,
                int,
                ENCODED_ELEM_T,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    this->alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    }
};

RCPP_MODULE(integer_kmers_computer) {
    Rcpp::class_<IntegerKMersComputer>("IntegerKMersComputer")
            .constructor<Rcpp::IntegerVector>()
            .method("addKMers", &IntegerKMersComputer::addKMers)
            .method("addGappedKMers", &IntegerKMersComputer::addGappedKMers)
            .method("toList", &IntegerKMersComputer::toList);

}
