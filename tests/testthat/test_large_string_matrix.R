library(testthat)
source("utils.R")

invoke_test <- function(expectedRes, ...) {
  res <- seqR::count_kmers(kmer_dictionary_name="unordered_map",
                           batch_size=100,
                           hash_dim=2,
                           verbose=FALSE,
                           parallel_mode=TRUE,
                           ...)
  expect_matrices_equal(as.matrix(res), expectedRes)
}

nrow <- 100
ncol <- 1000000
sequenceMatrix <-
  matrix(rep(rep("A", ncol), nrow), byrow = TRUE, nrow = nrow)

test_that("large string matrix with one elements", {
  expectedRes <- matrix(rep(999996, nrow),
                        nrow = nrow)
  colnames(expectedRes) <- c("A.A.A.A.A_0.0.0.0")
  
  invoke_test(
    expectedRes = expectedRes,
    alphabet = c("A"),
    sequences = sequenceMatrix,
    k = 5,
    positional = FALSE,
    with_kmer_counts = TRUE
  )
})
