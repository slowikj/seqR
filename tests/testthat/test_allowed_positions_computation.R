library(testthat)

invoke_test <- function(encoded_sequence, expected_result, print_res=FALSE) {
  res <- seqR::get_not_allowed_sequence_positions(encoded_sequence)
  if(print_res) {
    print(res)
  }
  expect_equal(res, expected_result)
}

test_that("an empty sequence gives a vector that contains only sentinels", {
  invoke_test(encoded_sequence=integer(0),
              expected_result=c(0,1))
})

test_that("a sequence with all allowed characters gives a vector that contains only sentinels", {
  invoke_test(encoded_sequence=c(1,2,3,4,5),
              expected_result=c(0,6))
})

test_that("a sequence (1,2,-1,3,4) gives c(0, 3, 6)", {
  invoke_test(encoded_sequence=c(1,2,-1,4,5),
              expected_result=c(0,3,6))
})

test_that("a sequence (1,2,-1,3,4,-1) gives c(0, 3, 6, 7)", {
  invoke_test(encoded_sequence=c(1,2,-1,3,4,-1),
              expected_result=c(0,3,6,7))
})
