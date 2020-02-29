library(testthat)

test_that("generates proper integer alphabet encoding", {
  input <- c(1,2,3)
  res <- seqR::encode_alphabet_integer(input)
  expect_setequal(names(res), input)
  expect_equal(0, sum(duplicated(res)))
})