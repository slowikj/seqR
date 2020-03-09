#ifndef KMER_STRINGS_CREATOR_H
#define KMER_STRINGS_CREATOR_H

#include<string>
#include<memory>
#include<Rcpp.h>
#include "dictionary.h"
#include "kmer_counter.h"
#include "hash/complex_hasher.h"

class KMerCreator {
public:
  KMerCreator(const Rcpp::StringVector& sequence,
              const Rcpp::IntegerVector& gaps,
              const std::string& itemSeparator):
  sequence(sequence),
  itemSeparator(itemSeparator) {
    this->gapsAccumulated = static_cast<Rcpp::IntegerVector>(Rcpp::cumsum(gaps + 1));
  }
  
  std::string get(int begin) const;
  
  std::string getPositional(int begin,
                            std::string positionSeparator) const;
  
private:
  const Rcpp::StringVector& sequence;
  std::string itemSeparator;
  Rcpp::IntegerVector gapsAccumulated;
  
  std::size_t getTotalSize(int begin, int separatorLength) const;
  
};

Rcpp::StringVector getKMerNames(
    const Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher>& kmerCountsDictionary,
    const Rcpp::StringVector& sequence,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer);

#endif
