#' @include validators.R
#' @include seqR_result.R

.count_contiguous_kmers_string_vector_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_contiguous_kmers_string_vector,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.count_contiguous_kmers_string_list_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_contiguous_kmers_string_list,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.count_gapped_kmers_string_vector_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_gapped_kmers_string_vector,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.count_gapped_kmers_string_list_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_gapped_kmers_string_list,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.invoke_kmer_function <- function(rcpp_counting_function, alphabet_validator, sequences, params) {
  alphabet <- params[["alphabet"]]
  alphabet_validator(alphabet)
  result_list <- rcpp_counting_function(sequences, alphabet, params)
  .prepare_final_result(result_list)
}

.prepare_final_result <- function(rcpp_result_list) {
  .convert_seqR_list_to_custom_matrix(rcpp_result_list)
}
