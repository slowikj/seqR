expect_matrices_equal <- function(a, b) {
  if(nrow(a) == 0 || ncol(a) == 0) {
    testthat::expect_true(nrow(a) == nrow(b) && ncol(a) == ncol(b))
  } else {
    testthat::expect_setequal(colnames(a), colnames(b))
    
    ordered_b <- b[, colnames(a)]
    testthat::expect_equal(as.vector(a), as.vector(ordered_b))
  }
}

convert_seqR_list_to_matrix <- function(seqR_list) {
  if(length(seqR_list$i) == 0) {
    matrix(nrow = seqR_list$seqNum, ncol=0)
  } else {
    as.matrix(slam::simple_triplet_matrix(
      i = seqR_list$i,
      j = seqR_list$j,
      v = seqR_list$v,
      dimnames = list(NULL, seqR_list$names)
    ))
  }
}

to_matrix <- function(v) {
  res <- matrix(unname(v), byrow=TRUE, nrow=1)
  colnames(res) <- names(v)
  res
}
