library(testthat)

invoke_test <- function(expected_res, ...) {
  res <- seqR::compute_polynomial_multihash(...)
  expect_equal(res, expected_res)
}

test_that("test 2 hashes for one item without position", {
  invoke_test(expected_res=c(1,1),
              P=c(101, 7),
              M=rep(1e9 + 33, 2),
              items=c(1),
              begin=0,
              position=-1)
})

test_that("test 2 hashes for two items without position", {
  invoke_test(expected_res=c(103, 9),
              P=c(101, 7),
              M=rep(1e9 + 33, 2),
              items=c(1, 2),
              begin=0,
              position=-1)
})

test_that("test 3 hashes for 3 items without position", {
  invoke_test(expected_res=c(10406, 66, 578),
              P=c(101, 7, 23),
              M=rep(1e9 + 33, 3),
              items=c(1,2,3),
              begin=0,
              position=-1)
})

test_that("test 3 hashes for 3 items with position 1", {
  invoke_test(expected_res=c(10406, 66, 578, 0),
              P=c(101, 7, 23),
              M=rep(1e9 + 33, 3),
              items=c(1,2,3),
              begin=0,
              position=1)
})

test_that("test 3 hashes for 3 items without the first one without position", {
  invoke_test(expected_res=c(205, 17, 49),
              P=c(101, 7, 23),
              M=rep(1e9 + 33, 3),
              items=c(1,2,3),
              begin=1,
              position=-1)
})

test_that("test 3 hashes for 3 items without the first one with position 2", {
  invoke_test(expected_res=c(205, 17, 49, 1),
              P=c(101, 7, 23),
              M=rep(1e9 + 33, 3),
              items=c(1,2,3),
              begin=1,
              position=2)
})