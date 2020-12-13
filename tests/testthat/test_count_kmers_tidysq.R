library(testthat)
source("utils.R")

invoke_test <- function(expected_res, ...) {
  res <- seqR::count_kmers(...)
  
  expect_matrices_equal(as.matrix(res), expected_res)
}

test_that("(tidysq input) count 2-mers", {
  sequences <- tidysq::as.sq(c("AAAAAC", "AAA", "AAAC"))
  expected_res <- matrix(c(
    4, 1,
    2, 0,
    2, 1), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A_0", "A.C_0")
  invoke_test(expected_res = expected_res,
              k = 2,
              alphabet=c("A", "C"),
              sequences = sequences,
              positional = FALSE,
              with_kmer_counts=TRUE)
})

test_that("(tidysq input) count 3-mers some items are not allowed", {
  sequences <- tidysq::as.sq(c("AAAACAAAAC", "AACTAAAA", "AACTAAAAC"))
  expected_res <- matrix(c(
    4, 0,
    2, 1,
    2, 1
  ), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.0", "T.A.A_0.0")
  invoke_test(expected_res = expected_res,
              alphabet=c("A", "T"),
              sequences = sequences,
              k=3,
              positional = FALSE,
              with_kmer_counts=TRUE)
})

test_that("(tidysq input) the k-mer is longer than the sequence", {
  sequences <- tidysq::as.sq(c("AAAACAAAAC", "AACTAAAA", "AACTAAAAC"))
  expected_res <- matrix(nrow=3, ncol=0)
  invoke_test(expected_res = expected_res,
              alphabet=c("A", "T"),
              sequences = sequences,
              k = 100,
              positional = FALSE,
              with_kmer_counts=TRUE)
})