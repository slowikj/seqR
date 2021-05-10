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
  Rcpp::List dimnames = res["dimnames"];
  return {
      res["i"],
      res["j"],
      res["v"],
      dimnames[1],
      res["nrow"]};
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

// [[Rcpp::export(".merge_nonempty_kmer_results")]]
Rcpp::List mergeNonEmptyKMerResults(
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

  int leftRowsNum = leftRows.size();
  for (std::size_t i = 0; i < leftRowsNum; ++i) {
    rows[i] = leftRows[i];
    cols[i] = leftCols[i];
    counts[i] = leftCounts[i];
  }

  int rightRowsNum = rightRows.size();
  for (std::size_t i = 0; i < rightRowsNum; ++i) {
    rows[i + leftRowsNum] = rightRows[i] + leftRowsNum;
    cols[i + leftRowsNum] = rightKMerIndexMapper[rightCols[i] - 1];
    counts[i + leftRowsNum] = rightCounts[i];
  }

  for (const auto &mappedKMer : kMerMapper) {
    kMers[mappedKMer.second - 1] = mappedKMer.first;
  }

  auto res = Rcpp::List::create(
      Rcpp::Named("i") = rows,
      Rcpp::Named("j") = cols,
      Rcpp::Named("v") = counts,
      Rcpp::Named("nrow") = leftNRow + rightNRow,
      Rcpp::Named("ncol") = kMers.size(),
      Rcpp::Named("dimnames") = Rcpp::List::create(R_NilValue, kMers));
  res.attr("class") = "simple_triplet_matrix";
  return res;
}
