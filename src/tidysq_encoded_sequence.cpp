#include "tidysq_encoded_sequence.h"
// [[Rcpp::depends(tidysq)]]
#include <tidysq.h>
#include <algorithm>
#include <iterator>
#include <RcppParallel.h>
#include <memory>

TidysqEncodedSequenceProxy::TidysqEncodedSequenceProxy(unsigned char* encodedSequence):
  encodedSequence(encodedSequence),
  size_(computeSize(encodedSequence)) {
}

TidysqEncodedSequenceProxy::TidysqEncodedSequenceProxy() :
  encodedSequence(nullptr),
  size_(0) {
}

TidysqEncodedSequenceProxy& TidysqEncodedSequenceProxy::operator=(const TidysqEncodedSequenceProxy& other) {
  if(this != &other) {
    this->setFields(other);
  }
  return (*this);
}

TidysqEncodedSequenceProxy::TidysqEncodedSequenceProxy(const TidysqEncodedSequenceProxy& other) {
  this->setFields(other);
}

TidysqEncodedSequenceProxy& TidysqEncodedSequenceProxy::operator=(TidysqEncodedSequenceProxy&& other) noexcept {
  if(this != &other) {
    this->setFields(other);
  }
  return (*this);
}

TidysqEncodedSequenceProxy::TidysqEncodedSequenceProxy(TidysqEncodedSequenceProxy&& other) noexcept {
  this->setFields(other);
}

void TidysqEncodedSequenceProxy::setFields(const TidysqEncodedSequenceProxy& other) noexcept {
  this->size_ = other.size_;
  this->encodedSequence = other.encodedSequence;
}
  
unsigned char TidysqEncodedSequenceProxy::operator[](int index) const {
  return this->encodedSequence[index];
}

unsigned char& TidysqEncodedSequenceProxy::operator[](int index) {
  return this->encodedSequence[index];
}

std::size_t TidysqEncodedSequenceProxy::size() const {
  return this->size_;
}

std::size_t TidysqEncodedSequenceProxy::computeSize(unsigned char* encodedSequence) {
  std::size_t currentIndex = 0;
  while(encodedSequence[currentIndex] != '\0') {
    ++currentIndex;
  }
  return currentIndex;
}

class TidysqUnpackingWorker: public RcppParallel::Worker {
public:
  TidysqUnpackingWorker(Rcpp::List& sq):
  sq(sq),
  alphabetSize(tidysq::get_alph_size(sq.attr("alphabet"))) {
    this->res.resize(sq.size());
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    for(int sq_i = begin; sq_i < end; ++sq_i) {
      res[sq_i] = std::move(TidysqEncodedSequenceProxy(
        tidysq::unpack_sq_to_char(sq[sq_i], alphabetSize)
      ));
    }
  }
  
private:
  Rcpp::List& sq;
  std::size_t alphabetSize;
  
public:
  std::vector<TidysqEncodedSequenceProxy> res;
};

std::vector<TidysqEncodedSequenceProxy> getEncodedTidysqSequences(Rcpp::List& sq) {
  TidysqUnpackingWorker unpackingWorker(sq);
  RcppParallel::parallelFor(0, sq.size(), unpackingWorker);
  return std::move(unpackingWorker.res);
}
