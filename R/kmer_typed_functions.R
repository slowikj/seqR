#' @include validators.R

count_kmers_integer_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_kmers_integer,
                        alphabet_validator=validate_integer,
                        sequences=sequences,
                        ...)
}

count_kmers_string_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_kmers_string,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

count_kmers_numeric_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_kmers_numeric,
                        alphabet_validator=validate_numeric,
                        sequences=sequences,
                        ...)
}

count_kmers_tidysq_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_kmers_tidysq,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_integer_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_gapped_kmers_integer,
                        alphabet_validator=validate_integer,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_string_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_gapped_kmers_string,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_numeric_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_gapped_kmers_numeric,
                        alphabet_validator=validate_numeric,
                        sequences=sequences,
                        ...)
}

count_gapped_kmers_tidysq_proxy <- function(sequences, ...) {
  .invoke_kmer_function(rcpp_counting_function=find_gapped_kmers_tidysq,
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        ...)
}

.invoke_kmer_function <- function(rcpp_counting_function, alphabet_validator, sequences, ...) {
  params <- list(...)
  alphabet_validator(params[["alphabet"]])
  result_list <- rcpp_counting_function(sequences, ...)
  slam::simple_triplet_matrix(
    i = result_list$i,
    j = result_list$j,
    v = result_list$v,
    dimnames = list(NULL, result_list$names)
  )
}
