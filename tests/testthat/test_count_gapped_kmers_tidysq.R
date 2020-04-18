library(testthat)
source("utils.R")

invoke_test <- function(expected_res, ...) {
  res <- seqR::count_gapped_kmers_tidysq(...)
  expect_matrices_equal(res, expected_res)
}

test_that("count non positional k-mers (0, 1)", {
  sq <- tidysq::construct_sq(c("AAAAAC", "AAA", "AAAC"), type="ami")
  expected_res <- matrix(c(
    2, 1,
    0, 0,
    0, 1), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A", "A.A.C")
  invoke_test(expected_res = expected_res,
              alphabet=c("A", "C"),
              sq = sq,
              gaps = c(0,1),
              positionalKMers = FALSE)
})
