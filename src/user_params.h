#pragma once

#include <string>
#include <vector>

#include "Rcpp.h"

struct UserParams {
  int k;
  std::vector<int> gaps;
  bool positional;
  bool withKMerCounts;
  bool withKMerNames;
  const std::string kMerDictionaryName;
  int batchSize;
  int hashDim;
  bool verbose;
  bool parallelMode;

  static UserParams createForContiguous(Rcpp::Environment &params) {
    UserParams res(params);
    res.gaps.resize(res.k - 1, 0);
    return res;
  }

  static UserParams createForGapped(Rcpp::Environment &params) {
    UserParams res(params);
    res.gaps = Rcpp::as<std::vector<int>>(params.get("gaps"));
    return res;
  }

 private:
  UserParams(Rcpp::Environment &params)
      : k(Rcpp::as<int>(params.get("k"))),
        positional(Rcpp::as<bool>(params.get("positional"))),
        withKMerCounts(Rcpp::as<bool>(params.get("with_kmer_counts"))),
        withKMerNames(Rcpp::as<bool>(params.get("with_kmer_names"))),
        kMerDictionaryName(Rcpp::as<std::string>(params.get("kmer_dictionary_name"))),
        batchSize(Rcpp::as<int>(params.get("batch_size"))),
        hashDim(Rcpp::as<int>(params.get("hash_dim"))),
        verbose(Rcpp::as<bool>(params.get("verbose"))),
        parallelMode(Rcpp::as<bool>(params.get("parallel_mode"))) {
  }
};
