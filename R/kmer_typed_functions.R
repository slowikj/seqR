#' @include validators.R

count_kmers_integer_proxy <- function(sequences, alphabet, ...) {
  validate_integer(alphabet)
  kmer_computer <- IntegerKMersComputer$new(alphabet)
  kmer_computer$addKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_kmers_string_proxy <- function(sequences, alphabet, ...) {
  validate_string(alphabet)
  kmer_computer <- StringKMersComputer$new(alphabet)
  kmer_computer$addKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_kmers_numeric_proxy <- function(sequences, alphabet, ...) {
  validate_numeric(alphabet)
  kmer_computer <- NumericKMersComputer$new(alphabet)
  kmer_computer$addKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_kmers_tidysq_proxy <- function(sequences, alphabet, ...) {
  validate_string(alphabet)
  tidysq_alphabet <- attr(sequences, "alphabet")
  kmer_computer <- TidySqKMersComputer$new(tidysq_alphabet, alphabet)
  kmer_computer$addKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_gapped_kmers_integer_proxy <- function(sequences, alphabet, ...) {
  validate_integer(alphabet)
  kmer_computer <- IntegerKMersComputer$new(alphabet)
  kmer_computer$addGappedKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_gapped_kmers_string_proxy <- function(sequences, alphabet, ...) {
  validate_string(alphabet)
  kmer_computer <- StringKMersComputer$new(alphabet)
  kmer_computer$addGappedKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_gapped_kmers_numeric_proxy <- function(sequences, alphabet, ...) {
  validate_numeric(alphabet)
  kmer_computer <- NumericKMersComputer$new(alphabet)
  kmer_computer$addGappedKMers(sequences, ...)
  .prepare_final_result(kmer_computer)
}

count_gapped_kmers_tidysq_proxy <- function(sequences, alphabet, ...) {
  validate_string(alphabet)
  tidysq_alphabet <- attr(sequences, "alphabet")
  kmer_computer <- TidySqKMersComputer$new(tidysq_alphabet, alphabet)
  print(list(...))
  kmer_computer$addGappedKMers(sq=sequences, ...)
  .prepare_final_result(kmer_computer)
}

.prepare_final_result <- function(kmer_computer) {
  result_list <- kmer_computer$toList()
  slam::simple_triplet_matrix(i=result_list$i,
                              j=result_list$j,
                              v=result_list$v,
                              dimnames=list(NULL, result_list$names))
}

