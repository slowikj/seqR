#' @include validators.R
#' @include seqR_result.R

.count_contiguous_kmers_string_vector_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_contiguous_kmers_string_vector,
                        kmer_alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.count_contiguous_kmers_string_list_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_contiguous_kmers_string_list,
                        kmer_alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.count_gapped_kmers_string_vector_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_gapped_kmers_string_vector,
                        kmer_alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.count_gapped_kmers_string_list_proxy <- function(sequences, params) {
  .invoke_kmer_function(rcpp_counting_function=.cpp_count_gapped_kmers_string_list,
                        kmer_alphabet_validator=validate_string,
                        sequences=sequences,
                        params=params)
}

.invoke_kmer_function <- function(rcpp_counting_function, kmer_alphabet_validator, sequences, params) {
  kmer_alphabet <- params[["kmer_alphabet"]]
  kmer_alphabet_validator(kmer_alphabet)
  result_list <- rcpp_counting_function(sequences, kmer_alphabet, params)
  .prepare_final_result(result_list)
}

.prepare_final_result <- function(rcpp_result_list) {
  .convert_internal_result_to_seqR_Matrix_class(rcpp_result_list)
}
