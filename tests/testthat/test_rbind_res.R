library(testthat)
source("utils.R")
library(seqR)

prepare_empty_seqR_res <- function(nrow) {
  Matrix::Matrix(ncol=0, nrow=nrow, sparse = TRUE)
}

prepare_seqR_res <- function(i, j, v, nrow, ncol, dimnames) {
  Matrix::sparseMatrix(i=i, j=j, x=v, dims = c(nrow, ncol), dimnames=dimnames)
}

test_that("merge two empty results", {
  resL <- prepare_empty_seqR_res(nrow=3)
  resR <- prepare_empty_seqR_res(nrow=5)
  
  res <- rbind_columnwise(resL, resR)
  expect_matrices_equal(as.matrix(res), matrix(nrow=8, ncol=0))
})

test_that("merge four empty results", {
  a <- prepare_empty_seqR_res(nrow=3)
  b <- prepare_empty_seqR_res(nrow=5)
  
  res <- rbind_columnwise(a, a, b, a)
  expect_matrices_equal(as.matrix(res), matrix(nrow=14, ncol=0))
})

test_that("merge empty and non empty results", {
  resL <- prepare_empty_seqR_res(nrow=4)
  resR <- prepare_seqR_res(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )

  res <- rbind_columnwise(resL, resR)
  
  expected_res <- prepare_seqR_res(
    i = c(5,5,5,5,6),
    j = c(1,2,3,1,2),
    v = c(1,2,3,4,5),
    nrow = 7,
    ncol = 3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )

  expect_matrices_equal(as.matrix(res), as.matrix(expected_res))
})

test_that("merge one non-empty and several empty results", {
  non_empty_res <- prepare_seqR_res(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  empty_res <- prepare_empty_seqR_res(nrow=5)
  
  res <- rbind_columnwise(non_empty_res, empty_res, empty_res)
  
  expected_res <- prepare_seqR_res(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=13,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  expect_matrices_equal(as.matrix(res), as.matrix(expected_res))
})

test_that("merge non empty and empty results", {
  resL <- prepare_seqR_res(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  resR <- prepare_empty_seqR_res(nrow=4)
  
  res <- rbind_columnwise(resL, resR)
  
  expected_res <- prepare_seqR_res(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=7,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  expect_matrices_equal(as.matrix(res), as.matrix(expected_res))
})

test_that("merge two non empty results with no k-mer in common", {
  resL <- prepare_seqR_res(
    i=c(1,1,3,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  resR <- prepare_seqR_res(
    i=c(3,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("X.X_1", "X.A_1", "X.B.A.A_1.1.0"))
  )
  
  expected_res <- prepare_seqR_res(
    i = c(1,1,3,1,2, 6,4,4,4,5),
    j = c(1,2,3,1,2, 4,5,6,4,5),
    v = c(1,2,3,4,5, 1,2,3,4,5),
    nrow = 6,
    ncol = 6,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1", "X.X_1", "X.A_1", "X.B.A.A_1.1.0"))
  )
  
  res <- rbind_columnwise(resL, resR)
  
  expect_matrices_equal(as.matrix(res), as.matrix(expected_res))
})

test_that("merge two same results", {
  resL <- prepare_seqR_res(
    i=c(1,3,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  resR <- resL
  
  expected_res <- prepare_seqR_res(
    i = c(1,3,1,1,2, 4,6,4,4,5),
    j = c(1,2,3,1,2, 1,2,3,1,2),
    v = c(1,2,3,4,5, 1,2,3,4,5),
    nrow=6,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  res <- rbind_columnwise(resL, resR)
  
  expect_matrices_equal(as.matrix(res), as.matrix(expected_res))
})

test_that("merge two non empty results that have one k-mer in common", {
  resL <- prepare_seqR_res(
    i=c(1,3,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  resR <- prepare_seqR_res(
    i=c(3,1,1,1,2),
    j=c(1,2,3,1,4),
    v=c(1,2,3,4,5),
    nrow=3,
    ncol=4,
    dimnames = list(NULL, c("X.X_1", "X.A_1", "X.B.A.A_1.1.0", "A.A_1"))
  )
  
  expected_res <- prepare_seqR_res(
    i=c(1,3,1,1,2,6,4,4,4,5),
    j=c(1,2,3,1,2,4,5,6,4,1),
    v=c(1,2,3,4,5,1,2,3,4,5),
    nrow=6,
    ncol=6,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1", "X.X_1", "X.A_1", "X.B.A.A_1.1.0"))
  )
  
  res <- rbind_columnwise(resL, resR)
  
  expect_matrices_equal(as.matrix(res), as.matrix(expected_res))
})
