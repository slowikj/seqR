// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include<memory>
#include "concrete_encoders.h"
#include "kmer_counter.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <string>
#include "kmer_strings_creator.h"
#include "kmer_counter.h"
#include "kmer_hash_indexer.h"

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

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix count_kmers(Rcpp::StringVector& alphabet,
                                Rcpp::StringMatrix& sequenceMatrix,
                                int k,
                                bool positionalKMers) {
  Rcpp::IntegerVector gaps(k-1);
  
  auto alphabetEncoding = getEncoding(alphabet);
  auto kmerCountsManagers = std::move(
    parallelComputeKMerCounts<Rcpp::StringMatrix,
                              Rcpp::StringMatrix::Row,
                              Rcpp::String::StringProxy,
                              std::string,
                              ENCODED_ELEM_T>(
        k,
        positionalKMers,
        sequenceMatrix,
        alphabetEncoding
    )
  );
  auto [hashIndexer, uniqueKMers] = indexKMerHashes(kmerCountsManagers);
  Rcpp::StringVector uniqueKMerStrings = std::move(
    parallelComputeKMerStrings(
      uniqueKMers,
      sequenceMatrix,
      gaps,
      positionalKMers,
      default_item_separator,
      default_position_separator
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
