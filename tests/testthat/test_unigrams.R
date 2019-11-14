library(testthat)

test_that("count positional unigrams", {
  res <- count_unigrams(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c("a", "b", "c"),
    pos=TRUE)
  
  expected <- c("0_c"=1, "2_a"=1, "1_b"=2, "0_a"=1, "2_c"=1)
  expect_mapequal(res, expected)
})

test_that("count non positional unigrams", {
  res <- count_unigrams(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c("a", "b", "c"),
    pos=FALSE)
  
  expected <- c("a"=2, "b"=2, "c"=2)
  expect_mapequal(res, expected)
})

test_that("count non positional unigrams if some elements are not in alphabet", {
  res <- count_unigrams(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c("a", "b"),
    pos=FALSE)
  
  expected <- c("a"=2, "b"=2)
  expect_mapequal(res, expected)
})

test_that("count positional unigrams if some elements are not in alphabet", {
  res <- count_unigrams(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c("a", "b"),
    pos=TRUE)
  
  expected <- c("0_a"=1, "1_b"=2, "2_a"=1) 
  expect_mapequal(res, expected)
})

test_that("print non existing in sequence unigrams with count 0 if pos is FALSE", {
  res <- count_unigrams(
    matrix(data=c("a", "a", "b", "b", "a", "b"), nrow=2),
    c("a", "b", "c"),
    pos=FALSE)
 
  expected <- c("a"=3, "b"=3, "c"=0)
  expect_mapequal(res, expected)
})