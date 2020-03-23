expect_matrices_equal <- function(a, b) {
  testthat::expect_setequal(colnames(a), colnames(b))
  
  ordered_b <- b[, colnames(a)]
  testthat::expect_equal(as.vector(a), as.vector(ordered_b))
}

to_matrix <- function(v) {
  res <- matrix(unname(v), byrow=TRUE, nrow=1)
  colnames(res) <- names(v)
  res
}