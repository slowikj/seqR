#ifndef TIDYSQ_ENCODED_SEQUENCE_H
#define TIDYSQ_ENCODED_SEQUENCE_H

#include <vector>
#include <Rcpp.h>

class TidysqEncodedSequenceProxy {
public:
  TidysqEncodedSequenceProxy();
  
  TidysqEncodedSequenceProxy(unsigned char* encodedSequence);
  
  TidysqEncodedSequenceProxy(const TidysqEncodedSequenceProxy&);
  
  TidysqEncodedSequenceProxy(TidysqEncodedSequenceProxy&&) noexcept;
  
  TidysqEncodedSequenceProxy& operator=(const TidysqEncodedSequenceProxy&);
  
  TidysqEncodedSequenceProxy& operator=(TidysqEncodedSequenceProxy&&) noexcept;
  
  unsigned char operator[](int index) const;
  
  unsigned char& operator[](int index);
  
  std::size_t size() const;
  
private:
  unsigned char* encodedSequence;
  
  std::size_t size_;
  
  std::size_t computeSize(unsigned char* encodedSequence);
  
  void setFields(const TidysqEncodedSequenceProxy& other) noexcept;
};

std::vector<TidysqEncodedSequenceProxy> getEncodedTidysqSequences(Rcpp::List& sq);

#endif
