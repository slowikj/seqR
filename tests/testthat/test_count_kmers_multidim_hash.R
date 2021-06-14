library(testthat)
source("utils.R")

invoke_test <- function(hash_dim, expected_res, ...) {
  res <- seqR::count_kmers(hash_dim = hash_dim, ...)
  expect_matrices_equal(expected_res, res)
}

invoke_dim_test <- function(hash_dim) {
  sq <- c("AAAAAABABA", "AAAABBAABABABABAB", "AAAA")
  expected_res <- matrix(c(
    1, 2, 4, 1, 0, 0, 0,
    4, 4, 2, 2, 1, 1, 1,
    0, 0, 2, 0, 0, 0, 0
  ), nrow=3, byrow=TRUE)
  colnames(expected_res) <- c("B.A.B_0.0", "A.B.A_0.0", "A.A.A_0.0", "A.A.B_0.0", "B.B.A_0.0", "A.B.B_0.0", "B.A.A_0.0")
  invoke_test(hash_dim = hash_dim,
              expected_res = expected_res,
              sequences = sq,
              k = 3,
              kmer_alphabet = c("A", "B"))
}

test_that("test count_kmers for 1 dimentional hash", {
  invoke_dim_test(1)
})

test_that("test count_kmers for 2 dimentional hash", {
  invoke_dim_test(2)
})

test_that("test count_kmers for 3 dimentional hash", {
  invoke_dim_test(3)
})

test_that("test count_kmers for 4 dimentional hash", {
  invoke_dim_test(4)
})

test_that("test count_kmers for 5 dimentional hash", {
  invoke_dim_test(5)
})

test_that("test count_kmers for 6 dimentional hash", {
  invoke_dim_test(6)
})

test_that("test count_kmers for 7 dimentional hash", {
  invoke_dim_test(7)
})

test_that("test count_kmers for 8 dimentional hash", {
  invoke_dim_test(8)
})

test_that("test count_kmers for 8 dimentional hash", {
  invoke_dim_test(50)
})

test_that("expect error for the hash dimension larger than 500", {
  expect_error(invoke_dim_test(501),
               info="hash_dim is a single integer number from the range [1, 500]")
})
