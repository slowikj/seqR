#include "tidysq_encoded_sequence.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include <algorithm>
#include <iterator>

std::vector<Rcpp::RawVector> getEncodedTidysqSequences(Rcpp::List& sq) {;
  auto alphabetSize = tidysq::get_alph_size(sq.attr("alphabet"));
  std::vector<Rcpp::RawVector> res;
  res.reserve(sq.size());
  std::transform(std::begin(sq), std::end(sq), std::back_inserter(res),
                 [&alphabetSize](const Rcpp::RawVector& rawSequence) -> Rcpp::RawVector {
                   return std::move(tidysq::unpack_raws(rawSequence, alphabetSize));
                 });
  return res;
}
