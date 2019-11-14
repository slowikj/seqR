library(testthat)

# hash for plain vector
test_that("get hash for 1 gives 1", {
  expect_equal(seqR::get_hash_for_word(c(1)), 1)
})

test_that("get hash for 5 gives 5", {
  expect_equal(seqR::get_hash_for_word(5), 5)
})

test_that("get hash for 10 gives 233", {
  expect_equal(seqR::get_hash_for_word(c(1, 0)), 233)
})

test_that("get hash for 201 gives 233^2 * 2 + 233^1 * 0 + 233^0 * 1", {
  expect_equal(seqR::get_hash_for_word(c(2, 0, 1)), 108579)
})

# hash for subsequence of vector
test_that("get hash for subword (not positional) gives valid result", {
  expect_equal(seqR::get_hash(c(1, 2, 3, 2, 0, 1), d=c(0,0), 3, pos=FALSE), 108579)
})

test_that("get hash for subword (positional) gives valid result - takes into account position in hash", {
  expect_equal(seqR::get_hash(c(1, 2, 3, 1, 2, 3), d=c(0,0), 1, pos=TRUE), 12758615)
})

test_that("get hash for subword (positional) that starts at 0 gives valid result", {
  expect_equal(seqR::get_hash(c(1, 2, 3, 1, 2, 3), d=c(0,0), 0, pos=TRUE), 54758)
})