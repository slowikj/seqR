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

class NumericKMersComputer {
public:
    using ENCODED_ELEM_T = short;

    explicit NumericKMersComputer(Rcpp::NumericVector alphabet) {
        initAlphabetEncoding(alphabet);
    }

    void addKMers(Rcpp::NumericMatrix sequenceMatrix,
                  int k,
                  bool positionalKMers,
                  bool withKMerCounts,
                  const std::string kmerDictionaryName,
                  int batchSize) {
        std::vector<int> gaps(k - 1);
        std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
        addKMersCommon(sequenceMatrix, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, algorithmParams,
                       batchSize,
                       result);
    }

    void addGappedKMers(Rcpp::NumericMatrix sequenceMatrix,
                        Rcpp::IntegerVector gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string kmerDictionaryName,
                        int batchSize) {
        auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
        auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
        addKMersCommon(sequenceMatrix, gapsConverted, positionalKMers, withKMerCounts, kmerDictionaryName,
                       hasherConfigs, batchSize, result);
    }

    Rcpp::List toList() const {
        return result.toRcppList();
    }

private:
    KMerCountingResult result;
    AlphabetEncoding<double, ENCODED_ELEM_T, UnorderedMapWrapper> alphabetEncoding;

    inline void initAlphabetEncoding(Rcpp::NumericVector &alphabet) {
        this->alphabetEncoding = std::move(
                prepareAlphabetEncodingFromRcpp<Rcpp::NumericVector, double, ENCODED_ELEM_T, UnorderedMapWrapper>(
                        alphabet));
    }

    template<class algorithm_params_t>
    inline
    void addKMersCommon(Rcpp::NumericMatrix &sequenceMatrix,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        int batchSize,
                        KMerCountingResult &kMerCountingResult) {
        for (int seqBegin = 0; seqBegin < sequenceMatrix.nrow(); seqBegin += batchSize) {
            int seqEnd = std::min(seqBegin + batchSize, static_cast<int>(sequenceMatrix.nrow()));
            addKMersCommon<algorithm_params_t>(sequenceMatrix, seqBegin, seqEnd,
                                               gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
                                               algorithmParams, kMerCountingResult);
        }
    }

    template<class algorithm_params_t>
    inline
    void addKMersCommon(Rcpp::NumericMatrix &sequenceMatrix,
                        int rowBegin,
                        int rowEnd,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        KMerCountingResult &kMerCountingResult) {
        KMerTaskConfig<RcppParallel::RMatrix<double>::Row, double> kMerTaskConfig(
                (rowEnd - rowBegin),
                getRMatrixRowGetter<Rcpp::NumericMatrix, double>(sequenceMatrix, rowBegin),
                gaps,
                positionalKMers,
                withKMerCounts,
                getIntToStringConverter(),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                RcppParallel::RMatrix<double>::Row,
                double,
                ENCODED_ELEM_T,
                UnorderedMapWrapper,
                algorithm_params_t>(kMerTaskConfig,
                                    this->alphabetEncoding,
                                    kmerDictionaryName,
                                    algorithmParams,
                                    kMerCountingResult);
    }
};

RCPP_MODULE(numeric_kmers_computer) {
    Rcpp::class_<NumericKMersComputer>("NumericKMersComputer")
            .constructor<Rcpp::NumericVector>()
            .method("addKMers", &NumericKMersComputer::addKMers)
            .method("addGappedKmers", &NumericKMersComputer::addGappedKMers)
            .method("toList", &NumericKMersComputer::toList);
}