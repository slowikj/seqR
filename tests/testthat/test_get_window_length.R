library(testthat)

test_that("get_window_length returns 3 for d=c(1)", {
  expect_equal(seqR::get_window_length(c(1)), 3)
})

test_that("get_window_length returns 2 for d=c(0)", {
  expect_equal(seqR::get_window_length(c(0)), 2)
})

test_that("get_window_length returns 6 for d=c(1,2)", {
  expect_equal(seqR::get_window_length(c(1,2)), 6)
})

test_that("get_window_length returns 10 for d=rep(0, 9)", {
  expect_equal(seqR::get_window_length(rep(0, 9)), 10)
})