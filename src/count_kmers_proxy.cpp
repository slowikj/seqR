// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
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

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector count_kmers(Rcpp::StringVector& alphabet,
                                Rcpp::StringVector& sequence,
                                int k,
                                bool isPositionalKMer) {
  auto alphabetEncoding = getEncoding(alphabet);
  auto kmerCountsManager = std::move(
    countKMers<Rcpp::StringVector, Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T>(
        k, sequence, alphabetEncoding, isPositionalKMer)
  );
  
  const Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher>& kmersDict = kmerCountsManager.getDictionary();  
  Rcpp::IntegerVector res(kmersDict.size());
  int resIndex = 0;
  for(const auto& kmerDictElem: kmersDict) {
    res[resIndex++] = kmerDictElem.second.cnt;
  }
  
  std::vector<KMerCountsManager> kmerCounts;
  kmerCounts.push_back(std::move(kmerCountsManager));
  auto [hashIndexer, uniqueKMers] = indexKMerHashes(kmerCounts);

  Rcpp::StringMatrix sequences(1, sequence.size());
  for(int i = 0; i < sequence.size(); ++i) {
    sequences[0, i] = sequence(i);
  }

  Rcpp::IntegerVector gaps(k-1);
  res.names() = parallelComputeKMerStrings(
    uniqueKMers,
    sequences,
    gaps,
    isPositionalKMer,
    ".",
    "_");
  
  return res;
}

// //' @export
// // [[Rcpp::export]]
// void count(Rcpp::StringVector alphabet,
//            Rcpp::StringMatrix sequenceMatrix,
//            int k,
//            bool isPositionalKMer) {
//   auto alphabetEncoding = getEncoding(alphabet);
//   auto res = std::move(
//     parallelComputeKMerCounts<Rcpp::StringMatrix, Rcpp::StringMatrix::Row, Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T>(
//       k,
//       isPositionalKMer,
//       sequenceMatrix,
//       alphabetEncoding
//   ));
//   for(int i = 0; i < res.size(); ++i) {
//     Rcpp::Rcout << "seq " << i << std::endl;
//     for(const auto& elem: (res[i]).getDictionary()) {
//       Rcpp::Rcout << "beg: " << elem.second.seqStartPosition << " cnt: " << elem.second.cnt << std::endl;
//     }
//   }
//   
// }
