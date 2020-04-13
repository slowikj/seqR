#include "sequence_getter.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>

// SequenceGetter_t<Rcpp::StringVector> getTidySqRowGetter(Rcpp::List& sq) {
//   int alphabetSize = tidysq::get_alph_size(sq.attr("alphabet"));
//   return [&sq, alphabetSize](int rowNum) -> Rcpp::StringVector {
//     return tidysq::unpack_chars(sq[rowNum], alphabetSize, NA_STRING);
//   };
// }
