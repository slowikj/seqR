#ifndef KMER_STRINGS_CREATOR_H
#define KMER_STRINGS_CREATOR_H

#include<string>
#include<memory>
// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include <functional>
#include <algorithm>
#include "dictionary.h"
#include "hash/complex_hasher.h"
#include "input_to_string_item_converter.h"
#include "sequence_getter.h"
#include "kmer_hash_indexer.h"

const std::string default_item_separator = ".";
const std::string default_section_separator = "_";

Rcpp::IntegerVector getGapsAccumulated(const Rcpp::IntegerVector& gaps);

template <class input_vector_t, class input_elem_t>
class KMerStringCreatorForSequence {
public:
  KMerStringCreatorForSequence(input_vector_t&& sequence,
                               const Rcpp::IntegerVector& gapsAccumulated,
                               const std::string& itemSeparator,
                               InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter):
    sequence(std::move(sequence)),
    gapsAccumulated(gapsAccumulated),
    itemSeparator(itemSeparator),
    inputToStringItemConverter(inputToStringItemConverter) {
  }
  
  std::string get(int begin) const {
    int totalSize = getTotalSize(begin, itemSeparator.size());
    std::string res;
    res.reserve(totalSize);
    res += inputToStringItemConverter(sequence[begin]);
    for(const int& accGap: this->gapsAccumulated) {
      res += itemSeparator + inputToStringItemConverter(sequence[begin + accGap]);
    }
    return res;
  }
  
  std::string getPositional(int begin,
                            std::string sectionSeparator) const {
    std::string withoutPositionString = this->get(begin);
    return std::to_string(begin + 1) + sectionSeparator + withoutPositionString;
  }
  
private:
  input_vector_t sequence;
  std::string itemSeparator;
  Rcpp::IntegerVector gapsAccumulated;
  InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter;
  
  std::size_t getTotalSize(int begin, int separatorLength) const {
    return std::accumulate(
        std::begin(this->gapsAccumulated),
        std::end(this->gapsAccumulated),
        inputToStringItemConverter(sequence[begin]).size(),
        [&toStrConverter=std::as_const(this->inputToStringItemConverter),
         &begin,
         &separatorLength,
         &seq=std::as_const(this->sequence)](int r, int accGap) {
          return r + toStrConverter(seq[begin + accGap]).size() + separatorLength;
        }
    );
  }
  
};

template <class input_vector_t, class input_elem_t>
class KMerStringsCreatorWorker: public RcppParallel::Worker {
public:
  Rcpp::StringVector outputKMerStrings;
  
  KMerStringsCreatorWorker(const std::vector<KMerPositionInfo>& indexedKMers,
                           int sequencesNum,
                           SequenceGetter_t<input_vector_t> sequenceGetter,
                           const Rcpp::IntegerVector& gaps,
                           bool isPositionalKMer,
                           InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter,
                           std::string itemSeparator,
                           std::string sectionSeparator):
    indexedKMers(indexedKMers),
    itemSeparator(itemSeparator),
    sectionSeparator(sectionSeparator),
    outputKMerStrings(std::move(Rcpp::StringVector(indexedKMers.size()))) {
    Rcpp::IntegerVector gapsAccumulated = std::move(getGapsAccumulated(gaps));
    prepareKMerStringsCreators(sequencesNum, gapsAccumulated, inputToStringItemConverter, sequenceGetter, itemSeparator);
    prepareCreateKMerFunc(isPositionalKMer);
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    for(int i = begin; i < end; ++i) {
      this->outputKMerStrings[i] = std::move(createKMerFunc(
        this->indexedKMers[i].seqNum,
        this->indexedKMers[i].position));
    }
  }
  
private:
  const std::vector<KMerPositionInfo>& indexedKMers;
  const std::string& itemSeparator;
  const std::string& sectionSeparator;
  std::vector<KMerStringCreatorForSequence<input_vector_t, input_elem_t>> kmerStringCreators;
  std::function<std::string(int, int)> createKMerFunc;
  
  void prepareKMerStringsCreators(int sequencesNum,
                                  const Rcpp::IntegerVector& gapsAccumulated,
                                  InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter,
                                  SequenceGetter_t<input_vector_t> sequenceGetter,
                                  std::string itemSeparator) {
    this->kmerStringCreators.reserve(sequencesNum);
    for(int i = 0; i < sequencesNum; ++i) {
      auto seq = std::move(sequenceGetter(i));
      this->kmerStringCreators.emplace_back(std::move(seq), gapsAccumulated, itemSeparator, inputToStringItemConverter);
    }
  }
  
  void prepareCreateKMerFunc(bool isPositionalKMer) {
    this->createKMerFunc = isPositionalKMer ?
      static_cast<std::function<std::string(int, int)>>(
          [this](int seqNum, int pos) {
            return kmerStringCreators[seqNum].getPositional(pos, sectionSeparator);
      }) :
      static_cast<std::function<std::string(int, int)>>(
          [this](int seqNum, int pos) {
          return kmerStringCreators[seqNum].get(pos);
      });
  }
  
};

template<class input_vector_t, class input_elem_t>
Rcpp::StringVector parallelComputeKMerStrings(
  const std::vector<KMerPositionInfo>& indexedKMers,
  InputToStringItemConverter_t<input_elem_t> inputToStringItemConverter,
  int sequencesNum,
  SequenceGetter_t<input_vector_t> sequenceGetter,
  const Rcpp::IntegerVector& gaps,
  bool isPositionalKMer,
  std::string itemSeparator = default_item_separator,
  std::string sectionSeparator = default_section_separator) {
  
  KMerStringsCreatorWorker<input_vector_t, input_elem_t> worker(
      indexedKMers,
      sequencesNum,
      sequenceGetter,
      gaps,
      isPositionalKMer,
      inputToStringItemConverter,
      itemSeparator,
      sectionSeparator);
  RcppParallel::parallelFor(0, indexedKMers.size(), worker);
  return std::move(worker.outputKMerStrings);
}

#endif
