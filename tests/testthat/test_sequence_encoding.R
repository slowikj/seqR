library(testthat)

test_that("enumerate string vector", {
  res <- seqR::enumerate_string_sequence(sequence=c("a", "a", "bb", "aa", "b"),
                                         alphabet=c("a", "aa"))
  expect_equal(res, c(1, 1, -1, 2, -1))
})

test_that("return a vector filled with -1 for an empty (null) alphabet", {
  res <- seqR::enumerate_string_sequence(sequence=c("a", "b", "c"),
                                         alphabet=c())
  
  expect_equal(res, rep(-1, 3))
})

test_that("return an empty vector if sequence is empty (null)", {
  res <- seqR::enumerate_string_sequence(sequence=c(),
                                         alphabet=c("a"))
  
  expect_equal(res, integer(0))
})