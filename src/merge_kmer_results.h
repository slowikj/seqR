#pragma once

#include <Rcpp.h>

#include <tuple>

#include "common_config.h"

namespace resultsMerging {
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
    Rcpp::IntegerVector>
initResultIntVectors(std::size_t entriesNum) {
  return {
      Rcpp::IntegerVector(entriesNum),  // i
      Rcpp::IntegerVector(entriesNum),  // j
      Rcpp::IntegerVector(entriesNum)   // v
  };
}

inline std::size_t getIntVectorLength(Rcpp::List res) {
  Rcpp::Nullable<Rcpp::IntegerVector> rows = res["i"];
  if (rows.isNull()) {
    return 0;
  } else {
    Rcpp::IntegerVector x = rows.get();
    return x.size();
  }
}

inline std::size_t computeResultIntVecLength(Rcpp::List resList) {
  std::size_t resultLength = 0;
  for (int i = 0; i < resList.size(); ++i) {
    Rcpp::List elem = resList[i];
    resultLength += getIntVectorLength(elem);
  }
  return resultLength;
}

inline Rcpp::List mergeKMerResults(Rcpp::List resList) {
  std::size_t processedSeqNum = 0;
  std::size_t itemsOffset = 0;
  auto [rows, cols, counts] = initResultIntVectors(computeResultIntVecLength(resList));
  dictionary::MartinusRobinHoodDictionary<Rcpp::String, uint32_t> kMerMapper;

  for (std::size_t resList_i = 0; resList_i < resList.size(); ++resList_i) {
    Rcpp::List currentRes = resList[resList_i];
    auto [currentRows, currentCols, currentCounts, currentKMers, currentNRow] = getParams(currentRes);

    std::vector<uint32_t> colMapper(currentKMers.size());
    for (int i = 0; i < currentKMers.size(); ++i) {
      auto newKMer = currentKMers[i];
      if (!kMerMapper.isPresent(newKMer)) {
        uint32_t newIndex = kMerMapper.size() + 1;
        kMerMapper[newKMer] = newIndex;
      }
      colMapper[i] = kMerMapper[newKMer];
    }

    for (int i = 0; i < currentRows.size(); ++i) {
      std::size_t shiftedIndex = i + itemsOffset;
      rows[shiftedIndex] = currentRows[i] + processedSeqNum;
      cols[shiftedIndex] = colMapper[currentCols[i] - 1];
      counts[shiftedIndex] = currentCounts[i];
    }

    processedSeqNum += currentNRow;
    itemsOffset += currentRows.size();
  }

  Rcpp::StringVector kMers(kMerMapper.size());
  for (const auto &mappedKMer : kMerMapper) {
    kMers[mappedKMer.second - 1] = mappedKMer.first;
  }

  return Rcpp::List::create(
      Rcpp::Named(config::PROXY_ROWS_NAME) = rows,
      Rcpp::Named(config::PROXY_COLUMNS_NAME) = cols,
      Rcpp::Named(config::PROXY_VALUES_NAME) = counts,
      Rcpp::Named(config::PROXY_NROW) = processedSeqNum,
      Rcpp::Named(config::PROXY_NCOL) = kMers.size(),
      Rcpp::Named(config::PROXY_COLUMN_NAMES_NAME) = kMers);
}
}  // namespace resultsMerging
