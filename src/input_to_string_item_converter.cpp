#include "input_to_string_item_converter.h"
#include <utility>
#include <iomanip>
#include <sstream>

InputToStringItemConverter_t<Rcpp::StringVector::stored_type> getRcppStringProxyToStringConverter() {
  return [](const Rcpp::StringVector::stored_type& elem) -> std::string {
    return Rcpp::as<std::string>(elem);
  };
}

InputToStringItemConverter_t<int> getIntToStringConverter() {
  return [](const int& elem) -> std::string {
    return std::to_string(elem);
  };
}

InputToStringItemConverter_t<double> getDoubleToStringConverter(int decimalPrecision) {
  return [decimalPrecision](const double& elem) -> std::string {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(decimalPrecision);
    stream << elem;
    return stream.str();
  };
}

InputToStringItemConverter_t<unsigned char> getEncodedTidySqItemToStringConverter(Rcpp::StringVector& decoder) {
  return [&decoder](const unsigned char& elem) -> std::string {
    return Rcpp::as<std::string>(decoder[elem - 1]);
  };
}
