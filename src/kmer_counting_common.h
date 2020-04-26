#ifndef KMER_COUNTING_COMMON_H
#define KMER_COUNTING_COMMON_H

#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include <vector>
#include "sequence_getter.h"
#include "kmer_counts_manager.h"
#include "kmer_hash_indexer.h"
#include "alphabet_encoder.h"
#include "kmer_strings_creator.h"
#include "input_to_string_item_converter.h"
#include <memory>
#include <functional>
#include "result_creator.h"

extern const std::string default_item_separator;
extern const std::string default_section_separator;

template<class input_vector_t>
using CountingKMersProc_t = std::function<KMerCountsManager(input_vector_t &)>;

template<class input_vector_t>
class KMerCounterWorker : public RcppParallel::Worker {
public:
    KMerCounterWorker(int rowsNum,
                      CountingKMersProc_t<input_vector_t> countingKMersProc,
                      SequenceGetter_t<input_vector_t> sequenceGetter) :
            countingKMersProc(countingKMersProc),
            sequenceGetter(sequenceGetter) {
        kmerCounts.resize(rowsNum);
    }

    inline void operator()(size_t begin, size_t end) override {
        for (int rowNum = begin; rowNum < end; ++rowNum) {
            auto row = std::move(sequenceGetter(rowNum));
            kmerCounts[rowNum] = std::move(countingKMersProc(row));
        }
    }

private:
    CountingKMersProc_t<input_vector_t> countingKMersProc;
    SequenceGetter_t<input_vector_t> sequenceGetter;

public:
    std::vector<KMerCountsManager> kmerCounts;
};

template<class input_vector_t, class input_elem_t, class encoded_elem_t, class alphabet_hasher_t>
inline
std::vector<KMerCountsManager> parallelComputeKMerCounts(
        int rowsNum,
        CountingKMersProc_t<input_vector_t> countingProc,
        SequenceGetter_t<input_vector_t> sequenceGetter) {
    KMerCounterWorker<input_vector_t> worker(rowsNum, countingProc, sequenceGetter);
    RcppParallel::parallelFor(0, rowsNum, worker);
    return std::move(worker.kmerCounts);
}

template<class input_vector_t, class input_elem_t, class encoded_elem_t, class alphabet_hasher_t>
using ParallelKMerCountingProc_t = std::function<std::vector<KMerCountsManager>(
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> &,
        SequenceGetter_t<input_vector_t>)>;

template<class input_vector_t, class input_elem_t, class encoded_elem_t, class alphabet_hasher_t>
inline
Rcpp::IntegerMatrix getKMerCountsMatrix(
        AlphabetEncoding<input_elem_t, encoded_elem_t, alphabet_hasher_t> &alphabetEncoding,
        int sequencesNum,
        SequenceGetter_t<input_vector_t> sequenceGetter,
        const Rcpp::IntegerVector &gaps,
        bool positionalKMers,
        ParallelKMerCountingProc_t<input_vector_t, input_elem_t, encoded_elem_t, alphabet_hasher_t> parallelKMerCountingProc,
        InputToStringItemConverter_t<input_elem_t> inputToStringConverter) {

    auto kmerCountsManagers = std::move(parallelKMerCountingProc(alphabetEncoding, sequenceGetter));

    auto[hashIndexer, uniqueKMers] = indexKMerHashes(kmerCountsManagers);
    Rcpp::StringVector uniqueKMerStrings = std::move(
            parallelComputeKMerStrings<input_vector_t, input_elem_t>(
                    uniqueKMers,
                    inputToStringConverter,
                    sequencesNum,
                    sequenceGetter,
                    gaps,
                    positionalKMers,
                    default_item_separator,
                    default_section_separator
            )
    );

    KMerMatrixCreatorWorker matrixCreatorWorker(
            sequencesNum,
            uniqueKMerStrings.size(),
            kmerCountsManagers,
            hashIndexer,
            uniqueKMerStrings
    );
    RcppParallel::parallelFor(0, sequencesNum, matrixCreatorWorker);

    return matrixCreatorWorker.outputKMerCounts;
}

#endif
