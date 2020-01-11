library(testthat)

invoke_sequence_test <- function(test_fun, sequence, alphabet, expected_res) {
  res <- test_fun(sequence=sequence, alphabet=alphabet)
  expect_equal(res, expected_res)
}

# string vectors
invoke_string_sequence_test <- function(sequence, alphabet, expected_res) {
  invoke_sequence_test(test_fun=seqR::enumerate_string_sequence,
                       sequence=sequence,
                       alphabet=alphabet,
                       expected_res)
}

test_that("enumerate string vector", {
  invoke_string_sequence_test(sequence=c("a", "a", "bb", "aa", "b"),
                              alphabet=c("a", "aa"),
                              expected_res=c(1, 1, -1, 2, -1))
})

test_that("return a vector filled with -1 for an empty (null) alphabet", {
  invoke_string_sequence_test(sequence=c("a", "b", "c"),
                              alphabet=c(),
                              expected_res=rep(-1, 3))
})

test_that("return an empty vector if sequence is empty (null)", {
  invoke_string_sequence_test(sequence=c(),
                              alphabet=c("a"),
                              expected_res=integer(0))
})


test_that("enumerate an integer vector, treated as a string vector", {
  invoke_string_sequence_test(sequence=c(1,2,3,1),
                              alphabet=c(1,2),
                              expected_res=c(1,2,-1,1))
})

# integer vectors

invoke_integer_sequence_test <- function(sequence, alphabet, expected_res) {
  invoke_sequence_test(test_fun=seqR::enumerate_integer_sequence,
                       sequence=sequence,
                       alphabet=alphabet,
                       expected_res)
}

test_that("enumerate an integer vector", {
  invoke_integer_sequence_test(sequence=c(1,2,3,1),
                               alphabet=c(1,2),
                               expected_res=c(1,2,-1,1))
})

# numeric vectors

invoke_numeric_sequence_test <- function(sequence, alphabet, expected_res) {
  invoke_sequence_test(test_fun=seqR::enumerate_numeric_sequence,
                       sequence=sequence,
                       alphabet=alphabet,
                       expected_res=expected_res)
}

test_that("enumerate a numeric vector", {
  invoke_numeric_sequence_test(sequence=c(1.5, 2.2, 1.5, 3.2, 2),
                               alphabet=c(1.5, 2),
                               expected_res=c(1,-1,1,-1,2))
})