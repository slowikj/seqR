library(testthat)

invoke_test <- function(expected_res, ...) {
  res <- seqR::compute_polynomial_hash(...)
  expect_equal(expected_res, res)
}

test_that("test hashing one element", {
  invoke_test(expected_res=1,
              P=101,
              M=1e9 + 33,
              items=c(1),
              begin=0,
              position=-1)
})

test_that("test hashing two elements", {
  invoke_test(expected_res=102,
              P=101,
              M=1e9 + 33,
              items=c(1, 1),
              begin=0,
              position=-1)
})

test_that("test hashing elements without the first element", {
  invoke_test(expected_res=1,
              P=101,
              M=1e9 + 33,
              items=c(1, 1),
              begin=1,
              position=-1)
})

test_that("test hashing two different elements", {
  invoke_test(expected_res=103,
              P=101,
              M=1e9 + 33,
              items=c(1, 2),
              begin=0,
              position=-1)
})

test_that("test hashing two different elements without the first element", {
  invoke_test(expected_res=2,
              P=101,
              M=1e9 + 33,
              items=c(1, 2),
              begin=1,
              position=-1)
})

test_that("test hashing three different elements without the first element", {
  invoke_test(expected_res=205,
              P=101,
              M=1e9 + 33,
              items=c(1, 2, 3),
              begin=1,
              position=-1)
})

test_that("test hashing with position", {
  invoke_test(expected_res=20506,
              P=101,
              M=1e9 + 33,
              items=c(2, 1),
              begin=0,
              position=3)
})
