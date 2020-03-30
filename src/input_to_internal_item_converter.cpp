#include "input_to_internal_item_converter.h"

InputToInternalItemConverter_t<int, int> getIntToIntConverter() {
  return [](const int& elem) -> int { return elem; };
}

InputToInternalItemConverter_t<double, double> getDoubleToDoubleConverter() {
  return [](const double& elem) -> double { return elem; };
}

InputToInternalItemConverter_t<Rcpp::String::StringProxy, std::string> getRcppStringToStringConverter() {
  return [](const Rcpp::String::StringProxy& elem) -> std::string { return Rcpp::as<std::string>(elem); };
}
