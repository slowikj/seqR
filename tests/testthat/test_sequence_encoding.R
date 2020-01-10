library(testthat)

test_that("enumerate string vector", {
  res <- seqR::enumerate_string_sequence(sequence=c("a", "a", "bb", "aa", "b"),
                                         alphabet=c("a", "aa"))
  expect_equal(res, c(1, 1, -1, 2, -1))
})
