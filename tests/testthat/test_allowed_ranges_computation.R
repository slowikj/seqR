library(testthat)

invoke_test <- function(encoded_sequence, window_length, expected_result, print_res=FALSE) {
  res <- seqR::get_valid_sequence_ranges(encoded_sequence, window_length)
  if(print_res) {
    print(res)
  }
  expect_equal(res, expected_result)
}

test_that("a sequence with all allowed characters and 3 window size gives one range corresponding to the whole one", {
  invoke_test(encoded_sequence=c(1,2,3,4,5),
              window_length=3,
              expected_result=list(begin=c(1), end=c(3)))
})

test_that("a sequence (1,2,-1,3,4) and 1 window size gives 2 ranges", {
  invoke_test(encoded_sequence= c(1,2,-1,4,5),
              window_length=1,
              expected_result=list(begin=c(1, 4), end=c(2, 5)))
})

test_that("a sequence (1,2,-1,3,4) and 3 window gives NO range", {
  invoke_test(encoded_sequence= c(1,2,-1,3,4),
              window_length=3,
              expected_result=list(begin=integer(0), end=integer(0)))
})

test_that("a sequence (1,2,2,-1,1,1,-1,-1,-1,1) and 3 window gives only [1,1] range", {
  invoke_test(encoded_sequence= c(1,2,2,-1,1,1,-1,-1,-1,1),
              window_length=3,
              expected_result=list(begin=c(1), end=c(1)))
})
