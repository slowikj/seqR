library(testthat)

test_that("count positional unigrams", {
  res <- count_kmers(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c(),
    c("a", "b", "c"),
    pos=TRUE)
  
  expected <- c("0_c"=1, "2_a"=1, "1_b"=2, "0_a"=1, "2_c"=1)
  expect_mapequal(res, expected)
})

test_that("count non positional unigrams", {
  res <- count_kmers(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c(),
    c("a", "b", "c"),
    pos=FALSE)
  
  expected <- c("a"=2, "b"=2, "c"=2)
  expect_mapequal(res, expected)
})