#ifndef INPUT_TO_INTERNAL_ITEM_CONVERTER_H
#define INPUT_TO_INTERNAL_ITEM_CONVERTER_H

#include <string>
#include <functional>
#include <Rcpp.h>

template <class input_elem_t, class internal_elem_t>
using InputToInternalItemConverter_t = std::function<internal_elem_t(const input_elem_t&)>;

InputToInternalItemConverter_t<int, int> getIntToIntConverter();

InputToInternalItemConverter_t<double, double> getDoubleToDoubleConverter();

InputToInternalItemConverter_t<Rcpp::String::StringProxy, std::string> getRcppStringToStringConverter();

#endif
