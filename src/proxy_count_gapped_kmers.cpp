// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include <memory>
#include <functional>
#include "kmer_task_solver.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_contiguous_intervals_matrix(const Rcpp::IntegerVector &gaps) {
    auto intervals = getContiguousIntervals(gaps);
    Rcpp::IntegerMatrix res(intervals.size(), 2);
    for (int i = 0; i < intervals.size(); ++i) {
        res(i, 0) = intervals[i].first + 1;
        res(i, 1) = intervals[i].second + 1;
    }
    return res;
}

inline std::vector<PolynomialSingleHasherConfig> getGappedKMerHasherConfigs() {
    std::vector<PolynomialSingleHasherConfig> res;
    res.emplace_back(101, 1e9 + 7);
    res.emplace_back(97, 1e9 + 33);
    return std::move(res);
}

template<class sequences_t,
        class alphabet_t>
inline
Rcpp::IntegerMatrix findKMers(sequences_t &sequences,
                              alphabet_t &alphabet,
                              Rcpp::IntegerVector &gaps,
                              bool positionalKMers,
                              bool withKMerCounts,
                              const std::string &kmerDictionaryName) {
    auto gapsConverted = std::move(convertRcppVector<int, Rcpp::IntegerVector>(gaps));
    auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
    return findKMers<decltype(hasherConfigs)>(
            sequences, alphabet, gapsConverted, positionalKMers, withKMerCounts, kmerDictionaryName, hasherConfigs);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix
find_gapped_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                         Rcpp::StringVector &alphabet,
                         Rcpp::IntegerVector &gaps,
                         bool positionalKMers,
                         bool withKMerCounts,
                         const std::string &kmerDictionaryName) {
    return findKMers(sequenceMatrix, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                              Rcpp::IntegerVector &alphabet,
                                              Rcpp::IntegerVector &gaps,
                                              bool positionalKMers,
                                              bool withKMerCounts,
                                              const std::string &kmerDictionaryName) {
    return findKMers(sequenceMatrix, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                              Rcpp::NumericVector &alphabet,
                                              Rcpp::IntegerVector &gaps,
                                              bool positionalKMers,
                                              bool withKMerCounts,
                                              const std::string &kmerDictionaryName) {
    return findKMers(sequenceMatrix, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName);
}


//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_gapped_kmers_tidysq(Rcpp::List &sq,
                                             Rcpp::StringVector &alphabet,
                                             Rcpp::IntegerVector &gaps,
                                             bool positionalKMers,
                                             bool withKMerCounts,
                                             const std::string &kmerDictionaryName) {
    return findKMers(sq, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName);
}
