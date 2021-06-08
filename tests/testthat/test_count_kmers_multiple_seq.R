library(testthat)
source("utils.R")

invoke_test <- function(expected_res, ...) {
  res <- seqR::count_kmers(hash_dim = 2,
                           verbose=FALSE,
                           ...)
  expect_matrices_equal(as.matrix(res), expected_res)
}

test_that("(string list) count non-positional 2-mers for 2 sequences", {
  sequences <- list(
    c("a", "a", "b", "a", "a", "b"),
    c("b", "a", "b", "b", "b", "a"))
  expected_res <- matrix(
    c(2, 2, 1, 0,
      0, 1, 2, 2),
    byrow=TRUE,
    nrow=2
  )
  colnames(expected_res) <- c("a.a_0", "a.b_0", "b.a_0", "b.b_0")
  invoke_test(expected_res=expected_res,
              kmer_alphabet=c("a", "b"),
              sequences=sequences,
              k=2,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string list) count positional 2-mers for 2 sequences", {
  sequences <- list(
    c("a", "a", "b", "a", "a", "b"),
    c("b", "a", "b", "b", "b", "a"))
  expected_res <- matrix(c(
    1, 0, 1, 1, 0, 1, 0, 1, 0,
    0, 1, 1, 0, 1, 0, 1, 0, 1
  ), byrow=TRUE, nrow=2)
  colnames(expected_res) <- c("1_a.a_0", "1_b.a_0", "2_a.b_0", "3_b.a_0", "3_b.b_0", "4_a.a_0", "4_b.b_0", "5_a.b_0", "5_b.a_0")
  invoke_test(expected_res=expected_res,
              kmer_alphabet=c("a", "b"),
              sequences=sequences,
              k=2,
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string list) count non positional 5-mers for 10 sequences (10^6 each)", {
  nrow <- 10
  ncol <- 1000000
  sequences <- lapply(1:nrow, function(i) rep("A", ncol))
  
  expected_res <- matrix(rep(999996, nrow),
                         nrow=nrow)
  colnames(expected_res) <- c("A.A.A.A.A_0.0.0.0")
  invoke_test(expected_res=expected_res,
              kmer_alphabet=c("A"),
              sequences=sequences,
              k=5,
              positional=FALSE,
              with_kmer_counts=TRUE)
})
