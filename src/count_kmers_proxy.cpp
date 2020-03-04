// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<memory>
#include "alphabet_encoder/concrete_encoders.h"
#include "kmer_counter.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector count_kmers(Rcpp::StringVector& alphabet,
                                Rcpp::StringVector& sequence,
                                int k,
                                bool isPositionalKMer) {
  auto alphabetEncoding = getEncoding(alphabet);
  auto kmerCountsManager = std::move(
    countKMers<Rcpp::StringVector, Rcpp::String::StringProxy, std::string, ENCODED_ELEM_T>(
        k, sequence, alphabetEncoding, isPositionalKMer)
  );
  return kmerCountsManager.getDictionary().size();
}