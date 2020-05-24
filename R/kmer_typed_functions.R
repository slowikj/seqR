#' @include validators.R

count_kmers_integer_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_kmers_integer(alphabet),
                        alphabet_validator=validate_integer,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_kmers_string_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_kmers_string(alphabet),
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_kmers_numeric_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_kmers_numeric(alphabet),
                        alphabet_validator=validate_numeric,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_kmers_tidysq_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_kmers_tidysq(alphabet, attr(sequences, "alphabet")),
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_gapped_kmers_integer_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_gapped_kmers_integer(alphabet),
                        alphabet_validator=validate_integer,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_gapped_kmers_string_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_gapped_kmers_string(alphabet),
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_gapped_kmers_numeric_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_gapped_kmers_numeric(alphabet),
                        alphabet_validator=validate_numeric,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

count_gapped_kmers_tidysq_proxy <- function(sequences, alphabet, ...) {
  .invoke_kmer_function(rcpp_counting_function=.prepare_find_gapped_kmers_tidysq(alphabet, attr(sequences, "alphabet")),
                        alphabet_validator=validate_string,
                        sequences=sequences,
                        alphabet=alphabet,
                        ...)
}

.prepare_find_kmers_integer <- function(sequences, alphabet) {
  .prepare_counting_kmer_function(IntegerKMersComputer, alphabet)
}

.prepare_find_gapped_kmers_integer <- function(sequences, alphabet) {
  .prepare_counting_gapped_kmer_function(IntegerKMersComputer, alphabet)
}

.prepare_find_kmers_string <- function(sequences, alphabet) {
  .prepare_counting_kmer_function(StringKMersComputer, alphabet)
}

.prepare_find_gapped_kmers_string <- function(sequences, alphabet) {
  .prepare_counting_gapped_kmer_function(StringKMersComputer, alphabet)
}

.prepare_find_kmers_numeric <- function(sequences, alphabet) {
  .prepare_counting_kmer_function(NumericKMersComputer, alphabet)
}

.prepare_find_gapped_kmers_numeric <- function(sequences, alphabet) {
  .prepare_counting_gapped_kmer_function(NumericKMersComputer, alphabet)
}

.prepare_find_kmers_tidysq <- function(sequences, alphabet, internal_alphabet_encoding) {
  .prepare_counting_kmer_function(sequences, alphabet, internal_alphabet_encoding)
}

.prepare_find_gapped_kmers_tidysq <- function(sequences, alphabet, internal_alphabet_encoding) {
  .prepare_counting_gapped_kmer_function(sequences, alphabet, internal_alphabet_encoding)
}

.prepare_counting_kmer_function <- function(init_object, ...) {
  computer <- .get_computer_object(init_object, ...)
  function(sequences, ...) {
    computer$addKMers(sequences, ...)
    computer
  }
}

.prepare_counting_gapped_kmer_function <- function(init_object, ...) {
  computer <- .get_computer_object(init_object, ...)
  function(sequences, ...) {
    computer$addGappedKMers(sequences, ...)
    computer
  }
}

.get_computer_object <- function(init_object, ...) {
  init_object$new(...)
}

.invoke_kmer_function <- function(rcpp_counting_function, alphabet_validator, sequences, alphabet, ...) {
  alphabet_validator(alphabet)
  r <- rcpp_counting_function(sequences, ...)
  .prepare_final_result(r)
}

.prepare_final_result <- function(r) {
  rcpp_result_list <- r$toList()
  if(length(rcpp_result_list$i) == 0) {
    slam::as.simple_triplet_matrix(matrix(nrow=rcpp_result_list$processed_sequences, ncol=0))
  } else {
    slam::simple_triplet_matrix(
      i = rcpp_result_list$i,
      j = rcpp_result_list$j,
      v = rcpp_result_list$v,
      dimnames = list(NULL, rcpp_result_list$names)
    )
  }
}
