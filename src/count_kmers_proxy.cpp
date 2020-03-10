// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<memory>
#include "concrete_encoders.h"
#include "kmer_counter.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include "kmer_strings_creator.h"

Rcpp::StringVector getKMerNames(
    const Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher>& kmerCountsDictionary,
    const Rcpp::StringVector& sequence,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer) {
  KMerCreator kmerCreator(sequence, gaps, ".");
  std::function<std::string(int)> createKMer = isPositionalKMer ?
  static_cast<std::function<std::string(int)>>([&kmerCreator](int begin) {
    return kmerCreator.getPositional(begin, "_");
  }) :
    static_cast<std::function<std::string(int)>>([&kmerCreator](int begin) {
      return kmerCreator.get(begin);
    });
  Rcpp::StringVector res(kmerCountsDictionary.size());
  int resIndex = 0;
  for(const auto& kmerDictElem: kmerCountsDictionary) {
    res[resIndex++] = createKMer(kmerDictElem.second.seqStartPosition);
  }
  return res;
}

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
  
  const Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher>& kmersDict = kmerCountsManager->getDictionary();  
  Rcpp::IntegerVector res(kmersDict.size());
  int resIndex = 0;
  for(const auto& kmerDictElem: kmersDict) {
    res[resIndex++] = kmerDictElem.second.cnt;
  }
  res.names() = getKMerNames(kmersDict, sequence, std::move(Rcpp::IntegerVector(k-1)), isPositionalKMer);
  
  return res;
}
