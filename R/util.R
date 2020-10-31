#' @export
convert_seqR_list_to_slam_matrix <- function(seqR_list) {
  if (length(seqR_list$i) == 0) {
    slam::as.simple_triplet_matrix(matrix(nrow = seqR_list$seqNum,
                                          ncol = 0))
  } else {
    slam::simple_triplet_matrix(
      i = seqR_list$i,
      j = seqR_list$j,
      v = seqR_list$v,
      nrow = seqR_list$seqNum,
      dimnames = list(NULL, seqR_list$names)
    )
  }
}