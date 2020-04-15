library(testthat)
source("utils.R")

invoke_test <- function(expected_res, ...) {
  res <- seqR::count_kmers_tidysq(...)
  expect_matrices_equal(res, expected_res)
}

test_that("count 3-mers for tidysq sequences A+", {
  sq <- tidysq::construct_sq(c("AAAAA", "AA", "AAAAAAA"), type = 'ami')
  expected_res <- matrix(c(
    3,
    0,
    5
  ), nrow=3)
  colnames(expected_res) <- c("A.A.A")
  invoke_test(expected_res=expected_res,
              alphabet=c("A"),
              sq=sq,
              k=3,
              positionalKMers=FALSE)
})

test_that("count 3-mers for tidysq sequences A+ longer", {
  sq <- tidysq::construct_sq(c(strrep("A", 1000000), strrep("A", 1000), strrep("A", 100)), type = 'ami')
  expected_res <- matrix(c(
    999998,
    998,
    98
  ), nrow=3)
  colnames(expected_res) <- c("A.A.A")
  invoke_test(expected_res=expected_res,
              alphabet=c("A"),
              sq=sq,
              k=3,
              positionalKMers=FALSE)
})