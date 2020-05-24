#include <Rcpp.h>
#include <string>
#include <algorithm>
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

    explicit IntegerKMersComputer(std::vector<int> alphabet) {
        initAlphabetEncoding(alphabet);
    }

    void addKMers(Rcpp::IntegerMatrix sequenceMatrix,
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

    void addGappedKMers(Rcpp::IntegerMatrix sequenceMatrix,
                        std::vector<int> gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string kmerDictionaryName,
                        int batchSize) {
        auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
        addKMersCommon(sequenceMatrix, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
                       hasherConfigs, batchSize, result);
    }

    Rcpp::List toList() const {
        return result.toRcppList();
    }

private:
    KMerCountingResult result;
    AlphabetEncoding<int, ENCODED_ELEM_T, UnorderedMapWrapper> alphabetEncoding;

    inline void initAlphabetEncoding(std::vector<int> &alphabet) {
        this->alphabetEncoding = std::move(
                getAlphabetEncoding<std::vector<int>, int, ENCODED_ELEM_T, UnorderedMapWrapper>(
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
    void addKMersCommon(Rcpp::IntegerMatrix &sequenceMatrix,
                        int rowBegin,
                        int rowEnd,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        KMerCountingResult &kMerCountingResult) {
        KMerTaskConfig<RcppParallel::RMatrix<int>::Row, int> kMerTaskConfig(
                (rowEnd - rowBegin),
                getRMatrixRowGetter<Rcpp::IntegerMatrix, int>(sequenceMatrix, rowBegin),
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
            .constructor<std::vector<int>>()
            .method("addKMers", &IntegerKMersComputer::addKMers)
            .method("addGappedKMers", &IntegerKMersComputer::addGappedKMers)
            .method("toList", &IntegerKMersComputer::toList);

}