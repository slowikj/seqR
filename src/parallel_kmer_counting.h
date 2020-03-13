// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include <vector>
#include "hash/complex_hasher.h"
#include "dictionary.h"
#include "kmer_counts_manager.h"
#include "kmer_counter.h"
#include <memory>
#include <functional>

template <class input_matrix_t, class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
class KMerCounterWorker : public RcppParallel::Worker {
public:
  std::vector<KMerCountsManager> kmerCounts;
  
  using CountingProcedure_t = std::function<KMerCountsManager(
    input_vector_t&,
    AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>&
  )>;

  KMerCounterWorker(input_matrix_t& sequenceMatrix,
                    AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding,
                    CountingProcedure_t countingKMersProc):
      sequenceMatrix(sequenceMatrix),
      alphabetEncoding(alphabetEncoding),
      countingKMersProc(countingKMersProc) {
    kmerCounts.resize(sequenceMatrix.nrow());
  }

  void operator()(size_t begin, size_t end) {
    for(int rowNum=begin; rowNum < end; ++rowNum) {
      auto row = this->sequenceMatrix(rowNum, Rcpp::_);
      kmerCounts[rowNum] = std::move(
        countingKMersProc(row, this->alphabetEncoding)
      );
    }
  }

private:
  input_matrix_t& sequenceMatrix;
  AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding;
  CountingProcedure_t countingKMersProc;

};

template <class input_matrix_t, class input_vector_t, class input_elem_t, class internal_elem_t, class encoded_elem_t>
std::vector<KMerCountsManager> parallelComputeKMerCounts(
    int k,
    bool positionalKMer,
    input_matrix_t& sequenceMatrix,
    AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& alphabetEncoding) {
  KMerCounterWorker<input_matrix_t, input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t> worker(
      sequenceMatrix,
      alphabetEncoding,
      [k, positionalKMer](
          input_vector_t& v,
          AlphabetEncoding<input_elem_t, internal_elem_t, encoded_elem_t>& enc
      ) -> KMerCountsManager {
        return countKMers<input_vector_t, input_elem_t, internal_elem_t, encoded_elem_t>(
            k,
            v,
            enc,
            positionalKMer
        );
      }
  );
  RcppParallel::parallelFor(0, sequenceMatrix.nrow(), worker);
  return std::move(worker.kmerCounts);
}
