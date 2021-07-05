#include "proxy_contiguous_kmers.h"
#include "proxy_gapped_kmers.h"
#include "proxy_merge_kmer_results.h"

//[[Rcpp::export(".cpp_count_contiguous_kmers_string_vector")]]
Rcpp::List count_contiguous_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &kmerAlphabet,
    Rcpp::Environment &rcppParams) {
  return countContiguousKMers(sq, kmerAlphabet, rcppParams);
}

// [[Rcpp::export(".cpp_count_contiguous_kmers_string_list")]]
Rcpp::List count_contiguous_kmers_string_list(
    Rcpp::List &sq,
    Rcpp::StringVector &kmerAlphabet,
    Rcpp::Environment &rcppParams) {
  return countContiguousKMers(sq, kmerAlphabet, rcppParams);
}

// [[Rcpp::export(".cpp_count_gapped_kmers_string_vector")]]
Rcpp::List count_gapped_kmers_string_vector(
    Rcpp::StringVector &sq,
    Rcpp::StringVector &kmerAlphabet,
    Rcpp::Environment &rcppParams) {
  return countGappedKMers(sq, kmerAlphabet, rcppParams);
}

// [[Rcpp::export(".cpp_count_gapped_kmers_string_list")]]
Rcpp::List count_gapped_kmers_string_list(
    Rcpp::List &sq,
    Rcpp::StringVector &kmerAlphabet,
    Rcpp::Environment &rcppParams) {
  return countGappedKMers(sq, kmerAlphabet, rcppParams);
}

// [[Rcpp::export(".cpp_merge_kmer_results")]]
Rcpp::List merge_kmer_results(Rcpp::List resList) {
  return resultsMerging::mergeKMerResults(resList);
}
