#' @include validators.R
#' @include util.R

count_kmers_integer_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_contiguous_kmers_integer,
                        alphabet_validator=validate_integer,
                        sequences=sequences,
                        ...)
}

count_kmers_string_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_contiguous_kmers_string,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

count_kmers_numeric_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_contiguous_kmers_numeric,
                        alphabet_validator=validate_numeric,
                        sequences=sequences,
                        ...)
}

count_kmers_list_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_contiguous_kmers_list,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_integer_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_gapped_kmers_integer,
                        alphabet_validator=validate_integer,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_string_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_gapped_kmers_string,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_numeric_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_gapped_kmers_numeric,
                        alphabet_validator=validate_numeric,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_list_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=count_gapped_kmers_list,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

.invoke_kmer_function <- function(rcpp_counting_function, alphabet_validator, sequences, ...) {
  params <- list(...)
  alphabet_validator(params[["alphabet"]])
  result_list <- rcpp_counting_function(sequences, ...)
  .prepare_final_result(result_list)
}

.prepare_final_result <- function(rcpp_result_list) {
  convert_seqR_list_to_slam_matrix(rcpp_result_list)
}
