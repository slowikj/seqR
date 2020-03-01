library(testthat)

is_integer_vector <- function(v) {
  all(v == as.integer(v))
}

validate_encoding <- function(encoding_fun, input) {
  res <- encoding_fun(input)
  expect_setequal(names(res), input)
  expect_equal(0, sum(duplicated(res)))
  expect_true(is_integer_vector(res))
}

test_that("generates proper integer alphabet encoding", {
  validate_encoding(encoding_fun=seqR::encode_integer_alphabet,
                    input=c(1,2,3))
})

test_that("generates proper numeric alphabet encoding", {
  validate_encoding(encoding_fun=seqR::encode_numeric_alphabet,
                    input=c(4.4, 3.3, 1.123))
})