expect_matrices_equal <- function(a, b) {
  testthat::expect_setequal(colnames(a), colnames(b))
  
  ordered_b <- b[, colnames(a)]
  testthat::expect_equal(as.vector(a), as.vector(ordered_b))
}
