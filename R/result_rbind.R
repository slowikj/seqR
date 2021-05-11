#' @include util.R
#' @export
rbind.seqR_simple_triplet_matrix <- function(...) {
  .convert_seqR_list_to_custom_matrix(.cpp_merge_kmer_results(list(...)))
}
