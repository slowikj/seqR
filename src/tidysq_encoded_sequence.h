#ifndef TIDYSQ_ENCODED_SEQUENCE_H
#define TIDYSQ_ENCODED_SEQUENCE_H

#include <vector>
#include <Rcpp.h>

class TidysqEncodedSequence {
public:
  TidysqEncodedSequence(unsigned char* encodedSequence);
  
  unsigned char operator[](int index) const;
  
  unsigned char& operator[](int index);
  
  std::size_t size() const;
  
private:
  unsigned char* encodedSequence;
  std::size_t size_;
  
  std::size_t computeSize(unsigned char* encodedSequence);
};

std::vector<TidysqEncodedSequence> getEncodedTidysqSequences(Rcpp::List& sq);

#endif
