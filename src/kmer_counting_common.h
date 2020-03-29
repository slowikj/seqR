#ifndef KMER_COUNTING_COMMON_H
#define KMER_COUNTING_COMMON_H

// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include <vector>
#include "kmer_counts_manager.h"
#include "kmer_hash_indexer.h"
#include "concrete_encoders.h"
#include "kmer_strings_creator.h"
#include <memory>
#include <functional>

extern const std::string default_item_separator;
extern const std::string default_section_separator;

template <class input_vector_t>
using CountingKMersProc_t = std::function<KMerCountsManager(input_vector_t&)>;

template <class input_vector_t>
using RowGetter_t = std::function<input_vector_t(int)>;

template <class input_vector_t>
class KMerCounterWorker : public RcppParallel::Worker {
public:
  KMerCounterWorker(int rowsNum,
                    CountingKMersProc_t<input_vector_t> countingKMersProc,
                    RowGetter_t<input_vector_t> rowGetter):
    countingKMersProc(countingKMersProc),
    rowGetter(rowGetter) {
    kmerCounts.resize(rowsNum);
  }
  
  void operator()(size_t begin, size_t end) {
    for(int rowNum=begin; rowNum < end; ++rowNum) {
      auto row = std::move(rowGetter(rowNum));
      kmerCounts[rowNum] = std::move(countingKMersProc(row));
    }
  }
  
private:
  CountingKMersProc_t<input_vector_t> countingKMersProc;
  RowGetter_t<input_vector_t> rowGetter;
  
public:
  std::vector<KMerCountsManager> kmerCounts;
};

template <class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<KMerCountsManager> parallelComputeKMerCounts(
    int rowsNum,
    CountingKMersProc_t<input_vector_t> countingProc,
    RowGetter_t<input_vector_t> rowGetter) {
  KMerCounterWorker<input_vector_t> worker(
      rowsNum,
      countingProc,
      rowGetter
  );
  RcppParallel::parallelFor(0, rowsNum, worker);
  return std::move(worker.kmerCounts);
}

class KMerMatrixCreatorWorker: public RcppParallel::Worker {
public:
  Rcpp::IntegerMatrix outputKMerCounts;
  
  KMerMatrixCreatorWorker(int nrow,
                          int ncol,
                          std::vector<KMerCountsManager>& kmerCountsManagers,
                          Dictionary<std::vector<int>, int, vector_int_hasher>& hashIndexer,
                          Rcpp::StringVector& uniqueKMerStrings):
    outputKMerCounts(std::move(Rcpp::IntegerMatrix(nrow, ncol))),
    kmerCountsManagers(kmerCountsManagers),
    hashIndexer(hashIndexer),
    uniqueKMerStrings(uniqueKMerStrings) {
    Rcpp::colnames(this->outputKMerCounts) = Rcpp::wrap(uniqueKMerStrings);
  }
  
  void operator()(std::size_t beginRow, std::size_t endRow) {
    for(int r = beginRow; r < endRow; ++r) {
      for(const auto& kmerHashPair: kmerCountsManagers[r].getDictionary()) {
        int c = hashIndexer[kmerHashPair.first];
        outputKMerCounts(r, c) = kmerHashPair.second.cnt;
      }
    }
  }
  
private:
  std::vector<KMerCountsManager>& kmerCountsManagers;
  Dictionary<std::vector<int>, int, vector_int_hasher>& hashIndexer;
  Rcpp::StringVector& uniqueKMerStrings;
  
};

template <class input_matrix_t>
Rcpp::IntegerMatrix getKMerCountsMatrix(
  input_matrix_t& sequenceMatrix,
  const Rcpp::IntegerVector& gaps,
  bool positionalKMers,
  std::function<std::vector<KMerCountsManager>()> parallelKMerCountingProc) {

  auto kmerCountsManagers = std::move(parallelKMerCountingProc());
  auto [hashIndexer, uniqueKMers] = indexKMerHashes(kmerCountsManagers);
  Rcpp::StringVector uniqueKMerStrings = std::move(
    parallelComputeKMerStrings(
      uniqueKMers,
      sequenceMatrix,
      gaps,
      positionalKMers,
      default_item_separator,
      default_section_separator
    )
  );
  KMerMatrixCreatorWorker worker(
      sequenceMatrix.nrow(),
      uniqueKMerStrings.size(),
      kmerCountsManagers,
      hashIndexer,
      uniqueKMerStrings
  );
  RcppParallel::parallelFor(0, sequenceMatrix.nrow(), worker);
  return worker.outputKMerCounts;
}

#endif
