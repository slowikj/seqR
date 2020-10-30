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
Rcpp::List findKMers(sequences_t &sequences,
                     int sequencesNum,
                     alphabet_t &alphabet,
                     std::vector<int> &gaps,
                     bool positionalKMers,
                     bool withKMerCounts,
                     const std::string &kmerDictionaryName,
                     int batchSize) {
    auto hasherConfigs = std::move(getGappedKMerHasherConfigs());
    return findKMers<sequences_t, alphabet_t, decltype(hasherConfigs)>(
            sequences, sequencesNum, alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName,
            hasherConfigs, batchSize);
}

//' @export
// [[Rcpp::export]]
Rcpp::List
find_gapped_kmers_string(Rcpp::StringMatrix &sequenceMatrix,
                         std::vector<std::string> &alphabet,
                         std::vector<int> &gaps,
                         bool positionalKMers,
                         bool withKMerCounts,
                         const std::string &kmerDictionaryName,
                         int batchSize) {
    return findKMers(sequenceMatrix, sequenceMatrix.nrow(), alphabet, gaps, positionalKMers, withKMerCounts,
                     kmerDictionaryName, batchSize);
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_gapped_kmers_integer(Rcpp::IntegerMatrix &sequenceMatrix,
                                     std::vector<int> &alphabet,
                                     std::vector<int> &gaps,
                                     bool positionalKMers,
                                     bool withKMerCounts,
                                     const std::string &kmerDictionaryName,
                                     int batchSize) {
    return findKMers(sequenceMatrix, sequenceMatrix.nrow(), alphabet, gaps, positionalKMers, withKMerCounts,
                     kmerDictionaryName, batchSize);
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_gapped_kmers_numeric(Rcpp::NumericMatrix &sequenceMatrix,
                                     std::vector<double> &alphabet,
                                     std::vector<int> &gaps,
                                     bool positionalKMers,
                                     bool withKMerCounts,
                                     const std::string &kmerDictionaryName,
                                     int batchSize) {
    return findKMers(sequenceMatrix, sequenceMatrix.nrow(), alphabet, gaps, positionalKMers, withKMerCounts,
                     kmerDictionaryName, batchSize);
}

//' @export
// [[Rcpp::export]]
Rcpp::List find_gapped_kmers_list(Rcpp::List &sq,
                                  Rcpp::StringVector &alphabet,
                                  std::vector<int> &gaps,
                                  bool positionalKMers,
                                  bool withKMerCounts,
                                  const std::string &kmerDictionaryName,
                                  int batchSize) {
    return findKMers(sq, sq.size(), alphabet, gaps, positionalKMers, withKMerCounts, kmerDictionaryName, batchSize);
}
