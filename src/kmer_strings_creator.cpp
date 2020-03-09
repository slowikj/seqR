#include "kmer_strings_creator.h"
#include<utility>
#include<functional>

std::size_t KMerCreator::getTotalSize(int begin, int separatorLength) const {
  return this->sequence[begin].size() + std::accumulate(
      std::begin(this->gapsAccumulated),
      std::end(this->gapsAccumulated),
      sequence[begin].size(),
      [&begin, &separatorLength, &seq=std::as_const(this->sequence)](int r, int accGap) {
        return r + seq[begin + accGap].size() + separatorLength;
      }
  );
}

std::string KMerCreator::get(int begin) const {
  int totalSize = getTotalSize(begin, itemSeparator.size());
  std::string res;
  res.reserve(totalSize);
  res += Rcpp::as<std::string>(sequence[begin]);
  for(const int& accGap: this->gapsAccumulated) {
    res += itemSeparator + Rcpp::as<std::string>(sequence[begin + accGap]);
  }
  return res;
}

std::string KMerCreator::getPositional(int position,
                                       std::string positionSeparator) const {
  std::string withoutPositionString = this->get(position);
  return std::to_string(position + 1) + positionSeparator + withoutPositionString;
}

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
