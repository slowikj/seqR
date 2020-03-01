library(testthat)

test_that("generates proper integer alphabet encoding", {
  input <- c(1,2,3)
  res <- seqR::encode_integer_alphabet(input)
  expect_setequal(names(res), input)
  expect_equal(0, sum(duplicated(res)))
})

# test_that("generates proper numeric alphabet encoding", {
#   input <- c("4.4", "3.3", "1.123")
#   res <- seqR::encode_alphabet_numeric(input)
#   expect_setequal(names(res), input)
#   expect_equal(0, sum(duplicated(res)))
#   expect_true(is.integer(res))
# })