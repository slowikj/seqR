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
  using CountingProcedure_t = std::function<KMerCountsManager(
    input_vector_t&
  )>;
  
  KMerCounterWorker(input_matrix_t& sequenceMatrix,
                    CountingProcedure_t countingKMersProc):
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
  CountingProcedure_t countingKMersProc;
  
public:
  std::vector<KMerCountsManager> kmerCounts;
};

#endif
