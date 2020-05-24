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
    explicit TidySqKMersComputer(std::vector<std::string> tidysqElementsEncoding,
                                 std::vector<std::string> alphabet) {
        initAlphabetEncoding(tidysqElementsEncoding, alphabet);
        this->safeElementsEncoding = tidysqElementsEncoding;
    }

    void addKMers(Rcpp::List sq,
                  int k,
                  bool positionalKMers,
                  bool withKMerCounts,
                  const std::string kmerDictionaryName,
                  int batchSize) {
        std::vector<int> gaps(k - 1);
        std::function<ComplexHasher()> algorithmParams = []() -> ComplexHasher { return createKMerComplexHasher(); };
        addKMersCommon(sq, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, algorithmParams, batchSize,
                       result);
    }

    void addGappedKMers(Rcpp::List sq,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string kmerDictionaryName,
                        int batchSize) {
        auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
        addKMersCommon(sq, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
                       hasherConfigs, batchSize, result);
    }

    Rcpp::List toList() const {
        return result.toRcppList();
    }

private:
    AlphabetEncoding<unsigned char, unsigned char, UnorderedMapWrapper> alphabetEncoding;
    std::vector<std::string> safeElementsEncoding;
    KMerCountingResult result;

    inline
    void initAlphabetEncoding(std::vector<std::string> &tidySqElementsEncoding,
                              std::vector<std::string> &alphabet) {
        this->alphabetEncoding = std::move(
                prepareAlphabetEncodingForTidysq<std::vector<std::string>, unsigned char, UnorderedMapWrapper>(
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
                        int batchSize,
                        KMerCountingResult &kMerCountingResult) {
        for (int seqBegin = 0; seqBegin < sq.size(); seqBegin += batchSize) {
            int seqEnd = std::min(seqBegin + batchSize, static_cast<int>(sq.size()));
            addKMersCommon(sq, seqBegin, seqEnd, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
                           algorithmParams, kMerCountingResult);
        }
    }

    template<class algorithm_params_t>
    inline
    void addKMersCommon(Rcpp::List &sq,
                        int seqBegin,
                        int seqEnd,
                        std::vector<int> &gaps,
                        bool positionalKMers,
                        bool withKMerCounts,
                        const std::string &kmerDictionaryName,
                        algorithm_params_t &algorithmParams,
                        KMerCountingResult &kMerCountingResult) {
        SafeSequencesWrapper<unsigned char> sequenceWrapper(std::move(getEncodedTidysqSequences(sq, seqBegin, seqEnd)));
        auto sequenceGetter = getSequenceGetter(sequenceWrapper);
        KMerTaskConfig<SafeSequencesWrapper<unsigned char>::Row, unsigned char> kMerTaskConfig(
                (seqEnd - seqBegin),
                sequenceGetter,
                gaps,
                positionalKMers,
                withKMerCounts,
                getEncodedTidySqItemToStringConverter(safeElementsEncoding),
                DEFAULT_KMER_ITEM_SEPARATOR,
                DEFAULT_KMER_SECTION_SEPARATOR);
        computeResult<
                SafeSequencesWrapper<unsigned char>::Row,
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
            .constructor<std::vector<std::string>, std::vector<std::string>>()
            .method("addKMers", &TidySqKMersComputer::addKMers)
            .method("addGappedKMers", &TidySqKMersComputer::addGappedKMers)
            .method("toList", &TidySqKMersComputer::toList);
}