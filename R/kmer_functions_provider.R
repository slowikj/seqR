.invoke_contiguous_kmer_function <- function(sequences, ...) {
  .get_kmer_function(sequences)[[1]](sequences=sequences, ...)
}

.invoke_gapped_kmer_function <- function(sequences, ...) {
  .get_kmer_function(sequences)[[2]](sequences=sequences, ...)
}

.get_kmer_function <- function(sequences) {
  if (is.vector(sequences) & is.character(sequences)) {
    .kmer_functions_map[["string_vector"]]
  } else if (is.list(sequences)) {
    .kmer_functions_map[["string_list"]]
  } else {
    stop("sequences param has unsupported type")
  }
}
  
#' @include kmer_typed_functions.R
.kmer_functions_map <- list(
  "string_list" = list(.count_contiguous_kmers_string_list_proxy, .count_gapped_kmers_string_list_proxy),
  "string_vector" = list(.count_contiguous_kmers_string_vector_proxy, .count_gapped_kmers_string_vector_proxy)
)
