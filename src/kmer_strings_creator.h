#ifndef KMER_STRINGS_CREATOR_H
#define KMER_STRINGS_CREATOR_H

#include<string>
#include<memory>
#include<Rcpp.h>
#include "dictionary.h"
#include "kmer_counter.h"
#include "hash/complex_hasher.h"

const std::string default_item_separator = ".";
const std::string default_position_separator = "_";

class KMerCreatorForSequence {
public:
  KMerCreatorForSequence(const Rcpp::StringVector& sequence,
                         const Rcpp::IntegerVector& gapsAccumulated,
                         const std::string& itemSeparator):
    sequence(sequence),
    gapsAccumulated(gapsAccumulated),
    itemSeparator(itemSeparator) {
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

Rcpp::IntegerVector getGapsAccumulated(const Rcpp::IntegerVector& gaps);

Rcpp::StringVector getKMerNames(
    const Dictionary<std::vector<int>, KMerHashInfo, vector_int_hasher>& kmerCountsDictionary,
    const Rcpp::StringVector& sequence,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator = default_item_separator,
    std::string positionSeparator = default_position_separator);

#endif
