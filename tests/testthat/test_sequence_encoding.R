library(testthat)

invoke_string_sequence_test <- function(sequence, alphabet, expected_res) {
  res <- seqR::enumerate_string_sequence(sequence=sequence, alphabet=alphabet)
  expect_equal(res, expected_res)
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
