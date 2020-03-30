#ifndef INPUT_ITEM_TO_STRING_CONVERTER_H
#define INPUT_ITEM_TO_STRING_CONVERTER_H

#include <string>
#include <Rcpp.h>
#include <functional>

template <class input_elem_t>
using InputToStringItemConverter_t = std::function<std::string(const input_elem_t&)>;

InputToStringItemConverter_t<Rcpp::String::StringProxy> getRcppStringProxyToStringConverter();

InputToStringItemConverter_t<int> getIntToStringConverter();

InputToStringItemConverter_t<double> getDoubleToStringConverter(int decimalPrecision);

#endif
