#ifndef SEQR_RCPP_TO_CPP_CONVERTERS_H
#define SEQR_RCPP_TO_CPP_CONVERTERS_H

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <memory>

template<class source_t, class target_t>
inline target_t convert(const source_t &source) {
    return Rcpp::as<target_t>(source);
}

template<>
inline int convert(const int &source) {
    return source;
}

template<>
inline double convert(const double &source) {
    return source;
}

template<class target_elem_t, class source_rcpp_vector_t>
inline std::vector<target_elem_t> convertRcppVector(const source_rcpp_vector_t &sourceVector) {
    std::vector<target_elem_t> res;
    std::transform(
            std::begin(sourceVector),
            std::end(sourceVector),
            std::back_inserter(res),
            [](const typename source_rcpp_vector_t::stored_type& sourceElem) -> target_elem_t {
                return convert<typename source_rcpp_vector_t::stored_type, target_elem_t>(sourceElem);
            });
    return std::move(res);
}

#endif //SEQR_RCPP_TO_CPP_CONVERTERS_H
