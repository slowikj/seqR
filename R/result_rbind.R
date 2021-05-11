#' @include util.R
#' @export
merge_kmer_results <- function(...) {
  .convert_seqR_list_to_slam_matrix(.merge_kmer_results(list(...)))
}