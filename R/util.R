.convert_seqR_list_to_slam_matrix <- function(seqR_list) {
  if (length(seqR_list$i) == 0) {
    slam::as.simple_triplet_matrix(matrix(nrow = seqR_list$seqNum,
                                          ncol = 0))
  } else {
    if(length(seqR_list$names) != 0) {
      dimnames <- list(NULL, seqR_list$names)
    } else {
      dimnames <- NULL
    }
    
    slam::simple_triplet_matrix(
      i = seqR_list$i,
      j = seqR_list$j,
      v = seqR_list$v,
      nrow = seqR_list$seqNum,
      dimnames = dimnames
    )
  }
}

.convert_seqR_list_to_custom_matrix <- function(seqR_list) {
  r <- .convert_seqR_list_to_slam_matrix(seqR_list)
  class(r) <- c("seqR_simple_triplet_matrix", class(r))
  r
}

is.seqR_simple_triplet_matrix <- function(x) {
  inherits(x, "seqR_simple_triplet_matrix")
}
