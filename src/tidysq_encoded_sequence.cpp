// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include "tidysq_encoded_sequence.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include <algorithm>
#include <iterator>
#include <vector>

std::vector<RcppParallel::RVector<unsigned char>> getEncodedTidysqSequences(Rcpp::List &sq) {
    auto alphabetSize = tidysq::C_get_alph_size(sq.attr("alphabet"));
    std::vector<RcppParallel::RVector<unsigned char>> res;
    for(int sq_i = 0; sq_i < sq.size(); ++sq_i) {
        Rcpp::RawVector packedSequence = sq[sq_i];
        Rcpp::RawVector unpackedSequence = tidysq::C_unpack_raws(packedSequence, alphabetSize);
        res.emplace_back(unpackedSequence);
    }
    return std::move(res);
}
