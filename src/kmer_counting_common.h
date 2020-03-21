#ifndef KMER_COUNTING_COMMON_H
#define KMER_COUNTING_COMMON_H

#include <RcppParallel.h>
#include <Rcpp.h>
#include <vector>
#include "kmer_counts_manager.h"
#include <memory>
#include <functional>

template <class input_matrix_t, class input_vector_t>
class KMerCounterWorker : public RcppParallel::Worker {
public:
  KMerCounterWorker(input_matrix_t& sequenceMatrix,
                    std::function<KMerCountsManager(input_vector_t&)> countingKMersProc):
    sequenceMatrix(sequenceMatrix),
    countingKMersProc(countingKMersProc) {
    kmerCounts.resize(sequenceMatrix.nrow());
  }
  
  void operator()(size_t begin, size_t end) {
    for(int rowNum=begin; rowNum < end; ++rowNum) {
      auto row = this->sequenceMatrix(rowNum, Rcpp::_);
      kmerCounts[rowNum] = std::move(countingKMersProc(row));
    }
  }
  
private:
  input_matrix_t& sequenceMatrix;
  std::function<KMerCountsManager(input_vector_t&)> countingKMersProc;
  
public:
  std::vector<KMerCountsManager> kmerCounts;
};

template <class input_matrix_t, class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<KMerCountsManager> parallelComputeKMerCounts(
    input_matrix_t& sequenceMatrix,
    std::function<KMerCountsManager(input_vector_t&)> countingProc) {
  KMerCounterWorker<input_matrix_t, input_vector_t> worker(
      sequenceMatrix,
      countingProc
  );
  RcppParallel::parallelFor(0, sequenceMatrix.nrow(), worker);
  return std::move(worker.kmerCounts);
}

#endif
