// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include "tidysq_encoded_sequence.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include <algorithm>
#include <iterator>

Rcpp::List getEncodedTidysqSequences(Rcpp::List &sq) {
    auto alphabetSize = tidysq::C_get_alph_size(sq.attr("alphabet"));
    Rcpp::List res(sq.size());
    for(int sq_i = 0; sq_i < sq.size(); ++sq_i) {
        Rcpp::RawVector packedSequence = sq[sq_i];
        res[sq_i] = tidysq::C_unpack_raws(packedSequence, alphabetSize);
    }
    return res;
}
