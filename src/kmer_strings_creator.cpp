#include "kmer_strings_creator.h"
#include <utility>
#include <iomanip>
#include <sstream>

Rcpp::IntegerVector getGapsAccumulated(const Rcpp::IntegerVector& gaps) {
  return static_cast<Rcpp::IntegerVector>(Rcpp::cumsum(gaps + 1));
}
