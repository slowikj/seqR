.invoke_contiguous_kmer_function <- function(sequences, ...) {
  .get_kmer_function(sequences)[[1]](sequences=sequences, ...)
}

.invoke_gapped_kmer_function <- function(sequences, ...) {
  .get_kmer_function(sequences)[[2]](sequences=sequences, ...)
}

.get_kmer_function <- function(sequences) {
  if (is.vector(sequences) & is.character(sequences)) {
    .kmer_functions_map[["string_vector"]]
  } else if (has_integers_only(sequences)) {
    .kmer_functions_map[["integer"]]
  } else if (is.numeric(sequences)) {
    .kmer_functions_map[["numeric"]]
  } else if (is.character(sequences)) {
    .kmer_functions_map[["string"]]
  } else {
    stop("sequences param has unsupported type")
  }
}
  
#' @include kmer_typed_functions.R
.kmer_functions_map <- list(
  "integer" = list(.count_contiguous_kmers_integer_proxy, .count_gapped_kmers_integer_proxy),
  "string" = list(.count_contiguous_kmers_string_proxy, .count_gapped_kmers_string_proxy),
  "numeric" = list(.count_contiguous_kmers_numeric_proxy, .count_gapped_kmers_numeric_proxy),
  "string_vector" = list(.count_contiguous_kmers_string_vector_proxy, .count_gapped_kmers_string_vector_proxy)
)
