#include "tidysq_encoded_sequence.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include <algorithm>
#include <iterator>

std::vector<TidysqEncodedSequence> getEncodedTidysqSequences(Rcpp::List& sq) {;
  auto alphabetSize = tidysq::get_alph_size(sq.attr("alphabet"));
  std::vector<TidysqEncodedSequence> res;
  res.reserve(sq.size());
  std::transform(std::begin(sq), std::end(sq), std::back_inserter(res),
                 [&alphabetSize](const Rcpp::RawVector& rawSequence) -> unsigned char* {
                   return tidysq::unpack_sq_to_char(rawSequence, alphabetSize);
                 });
  return res;
}

TidysqEncodedSequence::TidysqEncodedSequence(unsigned char* encodedSequence):
  encodedSequence(encodedSequence),
  size_(computeSize(encodedSequence)) {
}

unsigned char TidysqEncodedSequence::operator[](int index) const {
  return this->encodedSequence[index];
}

unsigned char& TidysqEncodedSequence::operator[](int index) {
  return this->encodedSequence[index];
}

std::size_t TidysqEncodedSequence::size() const {
  return this->size_;
}

std::size_t TidysqEncodedSequence::computeSize(unsigned char* encodedSequence) {
  std::size_t currentIndex = 0;
  while(encodedSequence[currentIndex] != '\0') {
    ++currentIndex;
  }
  return currentIndex;
}
