#ifndef TIDYSQ_ENCODED_SEQUENCE_H
#define TIDYSQ_ENCODED_SEQUENCE_H

#include <vector>
#include <Rcpp.h>

std::vector<std::vector<unsigned char>> getEncodedTidysqSequences(Rcpp::List &sq, int seqBegin, int seqEnd);

#endif
