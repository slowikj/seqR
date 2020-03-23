library(testthat)
source("utils.R")

invoke_test_string <- function(expectedRes, ...) {
  res <- seqR::count_kmers_string(...)
  expect_matrices_equal(res, expectedRes)
}

test_that("count non-positional 2-mers for 2 sequences", {
  sequenceMatrix <- matrix(
    c("a", "a", "b", "a", "a", "b",
      "b", "a", "b", "b", "b", "a"),
    byrow = TRUE,
    nrow = 2
  )
  expectedRes <- matrix(
    c(2, 2, 1, 0,
      0, 1, 2, 2),
    byrow=TRUE,
    nrow=2
  )
  colnames(expectedRes) <- c("a.a", "a.b", "b.a", "b.b")
  invoke_test_string(expectedRes=expectedRes,
                     alphabet=c("a", "b"),
                     sequenceMatrix=sequenceMatrix,
                     k=2,
                     positionalKMers=FALSE
  )
})
