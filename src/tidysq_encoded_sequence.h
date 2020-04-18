#ifndef TIDYSQ_ENCODED_SEQUENCE_H
#define TIDYSQ_ENCODED_SEQUENCE_H

#include <vector>
#include <Rcpp.h>

std::vector<Rcpp::RawVector> getEncodedTidysqSequences(Rcpp::List& sq);

#endif
