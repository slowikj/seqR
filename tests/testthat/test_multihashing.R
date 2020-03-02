library(testthat)

test_that("test 2 hashes for one item without position", {
  res <- seqR::compute_polynomial_multihash(P=c(101, 7),
                                            M=rep(1e9 + 33, 2),
                                            items=c(1),
                                            begin=0,
                                            position=-1)
  expect_equal(c(1, 1), res)
})

test_that("test 2 hashes for two items without position", {
  res <- seqR::compute_polynomial_multihash(P=c(101, 7),
                                            M=rep(1e9 + 33, 2),
                                            items=c(1, 2),
                                            begin=0,
                                            position=-1)
  expect_equal(c(103, 9), res)
})