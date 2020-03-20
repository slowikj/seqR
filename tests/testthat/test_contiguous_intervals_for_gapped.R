library(testthat)

test_that("test one contiguous interval", {
  expectedRes <- matrix(c(
    1, 5
  ), ncol=2)
  expect_equal(
    expectedRes,
    seqR::get_contiguous_intervals_matrix(gaps=rep(0, 4))
  )
})

test_that("test two contiguous intervals", {
  expectedRes <- matrix(c(
    1, 5,
    7, 8
  ), ncol=2, byrow=TRUE)
  expect_equal(
    expectedRes,
    seqR::get_contiguous_intervals_matrix(gaps=c(rep(0, 4), 1, 0))
  )
})

test_that("test all intervals of length 1", {
  expectedRes <- matrix(c(
    1, 1,
    3, 3,
    5, 5,
    7, 7
  ), ncol=2, byrow=TRUE)
  
  expect_equal(
    expectedRes,
    seqR::get_contiguous_intervals_matrix(gaps=c(1, 1, 1))
  )
})

test_that("test contiguous interval of length greater than 1 at the end", {
  expectedRes <- matrix(c(
    1, 1,
    3, 3,
    8, 11
  ), ncol=2, byrow=TRUE)
  
  expect_equal(
    expectedRes,
    seqR::get_contiguous_intervals_matrix(gaps=c(1, 4, 0, 0, 0))
  )
})

test_that("test two contiguous interval greater than 1 and singleton between", {
  expectedRes <- matrix(c(
    1, 5,
    7, 7,
    10, 15
  ), ncol=2, byrow=TRUE)
  
  expect_equal(
    expectedRes,
    seqR::get_contiguous_intervals_matrix(gaps=c(rep(0, 4), 1, 2, rep(0, 5)))
  )
})
