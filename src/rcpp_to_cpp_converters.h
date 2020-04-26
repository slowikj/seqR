#ifndef SEQR_RCPP_TO_CPP_CONVERTERS_H
#define SEQR_RCPP_TO_CPP_CONVERTERS_H

#include <Rcpp.h>

template<class source_t, class target_t>
inline target_t convert(const source_t& source) {
    return Rcpp::as<target_t>(source);
}

template<>
inline int convert(const int& source) {
    return source;
}

template<>
inline double convert(const double& source) {
    return source;
}

#endif //SEQR_RCPP_TO_CPP_CONVERTERS_H
