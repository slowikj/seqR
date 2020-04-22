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
  colnames(expectedRes) <- c("a.a_0", "a.b_0", "b.a_0", "b.b_0")
  invoke_test_string(expectedRes=expectedRes,
                     alphabet=c("a", "b"),
                     sequenceMatrix=sequenceMatrix,
                     k=2,
                     positionalKMers=FALSE
  )
})

test_that("count positional 2-mers for 2 sequences", {
  sequenceMatrix <- matrix(
    c("a", "a", "b", "a", "a", "b",
      "b", "a", "b", "b", "b", "a"),
    byrow = TRUE,
    nrow = 2
  )
  expectedRes <- matrix(c(
   1, 0, 1, 1, 0, 1, 0, 1, 0,
   0, 1, 1, 0, 1, 0, 1, 0, 1
  ), byrow=TRUE, nrow=2)
  colnames(expectedRes) <- c("1_a.a_0", "1_b.a_0", "2_a.b_0", "3_b.a_0", "3_b.b_0", "4_a.a_0", "4_b.b_0", "5_a.b_0", "5_b.a_0")
  invoke_test_string(expectedRes=expectedRes,
                     alphabet=c("a", "b"),
                     sequenceMatrix=sequenceMatrix,
                     k=2,
                     positionalKMers=TRUE)
})
