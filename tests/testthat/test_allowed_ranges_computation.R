library(testthat)

test_that("a sequence with all allowed characters and 3 window size gives one range corresponding to the whole one", {
  encoded_sequence <- c(1,2,3,4,5)
  window_length <- 3
  res <- seqR::get_valid_sequence_ranges(encoded_sequence, window_length)
  expect_equal(res,
               list(begin=c(1), end=c(3)))
})

test_that("a sequence (1,2,-1,3,4) and 1 window size gives 2 ranges", {
  encoded_sequence <- c(1,2,-1,4,5)
  window_length <- 1
  res <- seqR::get_valid_sequence_ranges(encoded_sequence, window_length)
  expect_equal(res,
               list(begin=c(1, 4), end=c(2, 5)))
})

test_that("a sequence (1,2,-1,3,4) and 3 window gives NO range", {
  encoded_sequence <- c(1,2,-1,3,4)
  window_length <- 3
  res <- seqR::get_valid_sequence_ranges(encoded_sequence, window_length)
  expect_equal(res,
               list(begin=integer(0), end=integer(0)))
})

test_that("a sequence (1,2,2,-1,1,1,-1,-1,-1,1) and 3 window gives only [1,1] range", {
  encoded_sequence <- c(1,2,2,-1,1,1,-1,-1,-1,1)
  window_length <- 3
  res <- seqR::get_valid_sequence_ranges(encoded_sequence, window_length)
  expect_equal(res,
                list(begin=c(1), end=c(1)))
})
