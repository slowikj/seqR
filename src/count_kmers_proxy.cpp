// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<memory>
#include "concrete_encoders.h"
#include "kmer_counter.h"
#include <iostream>
#include <algorithm>
#include "kmer_strings_creator.h"

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
  res.names() = getKMerNames(kmersDict, sequence, std::move(Rcpp::IntegerVector(k-1)), isPositionalKMer);
  
  return res;
}
