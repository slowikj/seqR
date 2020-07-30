// [[Rcpp::plugins("cpp17")]]
#include <Rcpp.h>
#include "tidysq_encoded_sequence.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include <algorithm>
#include <iterator>
#include <vector>

std::vector<std::vector<unsigned char>> getEncodedTidysqSequences(Rcpp::List &sq, int seqBegin, int seqEnd) {
    Rcpp::CharacterVector alphabet = sq.attr("alphabet");
    auto alphabetSize = tidysq::C_get_alph_size(alphabet);
    std::vector<std::vector<unsigned char>> res;
    for (int sq_i = seqBegin; sq_i < seqEnd; ++sq_i) {
        Rcpp::RawVector packedSequence = sq[sq_i];
        res.push_back(std::move(tidysq::unpack_raws_to_std_vector(
               packedSequence, alphabetSize)));
    }
    return std::move(res);
}
