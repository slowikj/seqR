#ifndef KMER_STRINGS_CREATOR_H
#define KMER_STRINGS_CREATOR_H

#include<string>
#include<memory>
// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
//' @importFrom  RcppParallel RcppParallelLibs
#include <RcppParallel.h>
#include<functional>
#include<algorithm>
#include "dictionary.h"
#include "hash/complex_hasher.h"
#include "kmer_hash_indexer.h"

const std::string default_item_separator = ".";
const std::string default_position_separator = "_";

Rcpp::IntegerVector getGapsAccumulated(const Rcpp::IntegerVector& gaps);

template <class input_vector_t, class input_elem_t>
class KMerStringCreatorForSequence {
public:
  KMerStringCreatorForSequence(input_vector_t&& sequence,
                               const Rcpp::IntegerVector& gapsAccumulated,
                               const std::string& itemSeparator,
                               std::function<std::string(const input_elem_t&)> inputItemToStringConverter):
    sequence(std::move(sequence)),
    gapsAccumulated(gapsAccumulated),
    itemSeparator(itemSeparator),
    inputItemToStringConverter(inputItemToStringConverter) {
  }
  
  std::string get(int begin) const {
    int totalSize = getTotalSize(begin, itemSeparator.size());
    std::string res;
    res.reserve(totalSize);
    res += inputItemToStringConverter(sequence[begin]);
    for(const int& accGap: this->gapsAccumulated) {
      res += itemSeparator + inputItemToStringConverter(sequence[begin + accGap]);
    }
    return res;
  }
  
  std::string getPositional(int begin,
                            std::string positionSeparator) const {
    std::string withoutPositionString = this->get(begin);
    return std::to_string(begin + 1) + positionSeparator + withoutPositionString;
  }
  
private:
  input_vector_t sequence;
  std::string itemSeparator;
  Rcpp::IntegerVector gapsAccumulated;
  std::function<std::string(const input_elem_t&)> inputItemToStringConverter;
  
  std::size_t getTotalSize(int begin, int separatorLength) const {
    return std::accumulate(
        std::begin(this->gapsAccumulated),
        std::end(this->gapsAccumulated),
        inputItemToStringConverter(sequence[begin]).size(),
        [&toStrConverter=std::as_const(this->inputItemToStringConverter),
         &begin,
         &separatorLength,
         &seq=std::as_const(this->sequence)](int r, int accGap) {
          return r + toStrConverter(seq[begin + accGap]).size() + separatorLength;
        }
    );
  }
  
};

template<class input_matrix_t, class input_vector_t, class input_elem_t>
class KMerStringsCreatorWorker: public RcppParallel::Worker {
public:
  Rcpp::StringVector outputKMerStrings;
  
  KMerStringsCreatorWorker(const std::vector<KMerPositionInfo>& indexedKMers,
                           input_matrix_t& sequences,
                           const Rcpp::IntegerVector& gaps,
                           bool isPositionalKMer,
                           std::function<std::string(const input_elem_t&)> inputItemToStringConverter,
                           std::string itemSeparator,
                           std::string positionSeparator):
    indexedKMers(indexedKMers),
    itemSeparator(itemSeparator),
    positionSeparator(positionSeparator),
    outputKMerStrings(std::move(Rcpp::StringVector(indexedKMers.size()))) {
    Rcpp::IntegerVector gapsAccumulated = std::move(getGapsAccumulated(gaps));
    prepareKMerStringsCreators(sequences, gapsAccumulated, inputItemToStringConverter, itemSeparator);
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
  const std::string& positionSeparator;
  std::vector<KMerStringCreatorForSequence<input_vector_t, input_elem_t>> kmerStringCreators;
  std::function<std::string(int, int)> createKMerFunc;
  
  void prepareKMerStringsCreators(input_matrix_t& sequences,
                                  const Rcpp::IntegerVector& gapsAccumulated,
                                  std::function<std::string(const input_elem_t&)> inputItemToStringConverter,
                                  std::string itemSeparator) {
    this->kmerStringCreators.reserve(sequences.nrow());
    for(int i = 0; i < sequences.nrow(); ++i) {
      auto seq = sequences(i, Rcpp::_);
      this->kmerStringCreators.emplace_back(std::move(seq), gapsAccumulated, itemSeparator, inputItemToStringConverter);
    }
  }
  
  void prepareCreateKMerFunc(bool isPositionalKMer) {
    this->createKMerFunc = isPositionalKMer ?
      static_cast<std::function<std::string(int, int)>>(
          [this](int seqNum, int pos) {
            return kmerStringCreators[seqNum].getPositional(pos, positionSeparator);
      }) :
      static_cast<std::function<std::string(int, int)>>(
          [this](int seqNum, int pos) {
          return kmerStringCreators[seqNum].get(pos);
      });
  }
  
};

template<class input_matrix_t, class input_vector_t, class input_elem_t>
Rcpp::StringVector parallelComputeKMerStrings(
  const std::vector<KMerPositionInfo>& indexedKMers,
  std::function<std::string(const input_elem_t&)> inputItemToStringConverter,
  input_matrix_t& sequences,
  const Rcpp::IntegerVector& gaps,
  bool isPositionalKMer,
  std::string itemSeparator = default_item_separator,
  std::string positionSeparator = default_position_separator) {
  
  KMerStringsCreatorWorker<input_matrix_t, input_vector_t, input_elem_t> worker(
      indexedKMers,
      sequences,
      gaps,
      isPositionalKMer,
      inputItemToStringConverter,
      itemSeparator,
      positionSeparator);
  RcppParallel::parallelFor(
    0,
    indexedKMers.size(),
    worker
  );
  return std::move(worker.outputKMerStrings);
}

Rcpp::StringVector parallelComputeKMerStrings(
    const std::vector<KMerPositionInfo>& indexedKMers,
    Rcpp::StringMatrix& sequences,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator = default_item_separator,
    std::string positionSeparator = default_position_separator);

Rcpp::StringVector parallelComputeKMerStrings(
    const std::vector<KMerPositionInfo>& indexedKMers,
    Rcpp::IntegerMatrix& sequences,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator = default_item_separator,
    std::string positionSeparator = default_position_separator);

Rcpp::StringVector parallelComputeKMerStrings(
    const std::vector<KMerPositionInfo>& indexedKMers,
    Rcpp::NumericMatrix& sequences,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator = default_item_separator,
    std::string positionSeparator = default_position_separator);

#endif
