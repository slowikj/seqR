#pragma once

namespace config {
using encoded_elem_t = unsigned char;

const std::string DEFAULT_KMER_ITEM_SEPARATOR = ".";

const std::string DEFAULT_KMER_SECTION_SEPARATOR = "_";

const std::string ALPHABET_ALL_LABEL = "all";

const Rcpp::String PROXY_ROWS_NAME = "i";
const Rcpp::String PROXY_COLUMNS_NAME = "j";
const Rcpp::String PROXY_VALUES_NAME = "v";
const Rcpp::String PROXY_COLUMN_NAMES_NAME = "names";
const Rcpp::String PROXY_NROW = "nrow";
const Rcpp::String PROXY_NCOL = "ncol";
}  // namespace config
