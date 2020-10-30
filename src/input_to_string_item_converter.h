#ifndef INPUT_ITEM_TO_STRING_CONVERTER_H
#define INPUT_ITEM_TO_STRING_CONVERTER_H

#include <string>
#include <Rcpp.h>
#include <functional>
#include <iomanip>
#include <sstream>

template<class input_elem_t>
using InputToStringItemConverter_t = std::function<std::string(const input_elem_t &)>;

inline InputToStringItemConverter_t<std::string> getStringToStringConverter() {
    return [](const std::string &str) -> std::string {
        return str;
    };
}

inline InputToStringItemConverter_t<char> getCharToStringConverter() {
    return [](const char &elem) -> std::string {
        return std::string(1, elem);
    };
}


inline InputToStringItemConverter_t<int> getIntToStringConverter() {
    return [](const int &elem) -> std::string {
        return std::to_string(elem);
    };
}

inline InputToStringItemConverter_t<double> getDoubleToStringConverter(int decimalPrecision) {
    return [decimalPrecision](const double &elem) -> std::string {
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(decimalPrecision);
        stream << elem;
        return stream.str();
    };
}


#endif
