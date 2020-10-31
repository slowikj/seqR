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
  as.matrix(seqR::convert_seqR_list_to_slam_matrix(seqR_list))
}

to_matrix <- function(v) {
  res <- matrix(unname(v), byrow=TRUE, nrow=1)
  colnames(res) <- names(v)
  res
}
