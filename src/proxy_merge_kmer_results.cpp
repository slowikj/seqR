// [[Rcpp::plugins("cpp17")]]

#include <Rcpp.h>

#include <tuple>

#include "dictionary/emilib_hash_map_wrapper.h"

inline std::tuple<
    Rcpp::IntegerVector,
    Rcpp::IntegerVector,
    Rcpp::IntegerVector,
    Rcpp::StringVector,
    std::size_t>
getParams(Rcpp::List res) {
  Rcpp::Nullable<Rcpp::List> dimnames = res["dimnames"];
  if (dimnames.isNull()) {
    return {
        Rcpp::IntegerVector(),
        Rcpp::IntegerVector(),
        Rcpp::IntegerVector(),
        Rcpp::StringVector(),
        res["nrow"]};
  } else {
    Rcpp::List nonNullDimnames = dimnames.get();
    return {
        res["i"],
        res["j"],
        res["v"],
        nonNullDimnames[1],
        res["nrow"]};
  }
}

inline std::tuple<
    Rcpp::IntegerVector,
    Rcpp::IntegerVector,
    Rcpp::IntegerVector,
    Rcpp::StringVector>
initResult(std::size_t entriesNum, std::size_t kMersNum) {
  return {
      Rcpp::IntegerVector(entriesNum),  // i
      Rcpp::IntegerVector(entriesNum),  // j
      Rcpp::IntegerVector(entriesNum),  // v
      Rcpp::StringVector(kMersNum)      // kMers
  };
}

// [[Rcpp::export(".merge_kmer_results")]]
Rcpp::List mergeKMerResults(
    Rcpp::List resLeft,
    Rcpp::List resRight) {
  auto [leftRows, leftCols, leftCounts, leftKMers, leftNRow] = getParams(resLeft);
  auto [rightRows, rightCols, rightCounts, rightKMers, rightNRow] = getParams(resRight);

  dictionary::EmilibHashMapWrapper<Rcpp::String, uint32_t> kMerMapper;
  for (std::size_t i = 0; i < leftKMers.size(); ++i) {
    kMerMapper[leftKMers[i]] = i + 1;
  }

  std::vector<uint32_t> rightKMerIndexMapper(rightKMers.size());
  for (std::size_t i = 0; i < rightKMers.size(); ++i) {
    if (!kMerMapper.isPresent(rightKMers[i])) {
      int newIndex = kMerMapper.size() + 1;
      kMerMapper[rightKMers[i]] = newIndex;
    }
    rightKMerIndexMapper[i] = kMerMapper[rightKMers[i]];
  }

  auto [rows, cols, counts, kMers] = initResult(
      leftRows.size() + rightRows.size(),
      kMerMapper.size());

  for (std::size_t i = 0; i < leftRows.size(); ++i) {
    rows[i] = leftRows[i];
    cols[i] = leftCols[i];
    counts[i] = leftCounts[i];
  }

  for (std::size_t i = 0; i < rightRows.size(); ++i) {
    int index = i + leftRows.size();
    rows[index] = rightRows[i] + leftNRow;
    cols[index] = rightKMerIndexMapper[rightCols[i] - 1];
    counts[index] = rightCounts[i];
  }

  for (const auto& mappedKMer : kMerMapper) {
    kMers[mappedKMer.second - 1] = mappedKMer.first;
  }

  return Rcpp::List::create(
      Rcpp::Named("i") = rows,
      Rcpp::Named("j") = cols,
      Rcpp::Named("v") = counts,
      Rcpp::Named("seqNum") = leftNRow + rightNRow,
      Rcpp::Named("names") = kMers);
}
