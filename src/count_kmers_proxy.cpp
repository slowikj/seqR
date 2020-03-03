// [[Rcpp::plugins("c++17")]]
#include<Rcpp.h>
#include<memory>
#include "alphabet_encoder/concrete_encoders.h"
//#include "kmer_counter.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector count_kmers(Rcpp::StringVector& alphabet,
                                Rcpp::StringVector& sequence,
                                int k,
                                bool isPositionalKMer) {
  // auto alphabetEncoding = getEncoding(alphabet);
  //auto kmerCountsManager = std::move(countKMers(k, sequence, alphabetEncoding, isPositionalKMer));
  //return kmerCountsManager.getDictionary().size();
  return 1;
}