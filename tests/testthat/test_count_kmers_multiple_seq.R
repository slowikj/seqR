library(testthat)

invoke_test <- function(expectedRes, ...) {
  res <- seqR::count_kmers(...)
  expect_setequal(colnames(expectedRes), colnames(res))
  expectedRes <- expectedRes[, colnames(res)]
  expect_equal(as.vector(expectedRes), as.vector(res))
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
  invoke_test(expectedRes=expectedRes,
              alphabet=c("a", "b"),
              sequenceMatrix=sequenceMatrix,
              k=2,
              positionalKMers=FALSE
  )
})
