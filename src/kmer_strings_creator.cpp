#include "kmer_strings_creator.h"
#include <utility>

Rcpp::IntegerVector getGapsAccumulated(const Rcpp::IntegerVector& gaps) {
  return static_cast<Rcpp::IntegerVector>(Rcpp::cumsum(gaps + 1));
}

Rcpp::StringVector parallelComputeKMerStrings(
    const std::vector<KMerPositionInfo>& indexedKMers,
    Rcpp::StringMatrix& sequences,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator,
    std::string positionSeparator) {
  return std::move(parallelComputeKMerStrings<Rcpp::StringMatrix,  Rcpp::StringMatrix::Row, Rcpp::String::StringProxy>(
      indexedKMers,
      [](const Rcpp::String::StringProxy& elem) -> std::string { return Rcpp::as<std::string>(elem); },
      sequences,
      gaps,
      isPositionalKMer,
      itemSeparator,
      positionSeparator
  )); 
}

Rcpp::StringVector parallelComputeKMerStrings(
    const std::vector<KMerPositionInfo>& indexedKMers,
    Rcpp::IntegerMatrix& sequences,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator,
    std::string positionSeparator) {
  return std::move(parallelComputeKMerStrings<Rcpp::IntegerMatrix,  Rcpp::IntegerMatrix::Row, int>(
      indexedKMers,
      [](const int& elem) -> std::string { return std::to_string(elem); },
      sequences,
      gaps,
      isPositionalKMer,
      itemSeparator,
      positionSeparator
  ));
}

Rcpp::StringVector parallelComputeKMerStrings(
    const std::vector<KMerPositionInfo>& indexedKMers,
    Rcpp::NumericMatrix& sequences,
    const Rcpp::IntegerVector& gaps,
    bool isPositionalKMer,
    std::string itemSeparator,
    std::string positionSeparator) {
  return std::move(parallelComputeKMerStrings<Rcpp::NumericMatrix,  Rcpp::NumericMatrix::Row, double>(
      indexedKMers,
      [](const double& elem) -> std::string { return std::to_string(elem); },
      sequences,
      gaps,
      isPositionalKMer,
      itemSeparator,
      positionSeparator
  ));
}
