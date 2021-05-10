#' @include util.R
#' @export
merge_kmer_results <- function(resL, resR) {
  .convert_seqR_list_to_slam_matrix(.merge_kmer_results(resL, resR))
}