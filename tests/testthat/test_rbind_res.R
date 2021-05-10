library(testthat)
source("utils.R")

test_that("merge two empty results", {
  resL <- slam::as.simple_triplet_matrix(matrix(ncol=0, nrow=3))
  resR <- slam::as.simple_triplet_matrix(matrix(ncol=0, nrow=5))
  
  res <- seqR::merge_kmer_results(resL, resR)
  expect_matrices_equal(as.matrix(res), matrix(nrow=8, ncol=0))
})

test_that("merge empty and non empty results", {
  resL <- slam::as.simple_triplet_matrix(matrix(ncol=0, nrow=4))
  resR <- slam::simple_triplet_matrix(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )

  res <- seqR::merge_kmer_results(resL, resR)
  
  expected_res <- resR
  expected_res$i <- expected_res$i + 4;
  expected_res$nrow <- 7

  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("merge non empty and empty results", {
  resL <- slam::simple_triplet_matrix(
    i=c(1,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  resR <- slam::as.simple_triplet_matrix(matrix(ncol=0, nrow=4))
  
  res <- seqR::merge_kmer_results(resL, resR)
  
  expected_res <- resL
  expected_res$nrow <- 7
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("merge two non empty results with no k-mer in common", {
  resL <- slam::simple_triplet_matrix(
    i=c(1,1,3,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  resR <- slam::simple_triplet_matrix(
    i=c(3,1,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    dimnames = list(NULL, c("X.X_1", "X.A_1", "X.B.A.A_1.1.0"))
  )
  
  expected_res <- slam::simple_triplet_matrix(
    i = c(resL$i, resR$i + resL$nrow),
    j = c(resL$j, resR$j + ncol(resL)),
    v = c(resL$v, resR$v),
    nrow = resL$nrow + resR$nrow,
    dimnames = list(NULL, c(resL$dimnames[[2]], resR$dimnames[[2]]))
  )
  
  res <- seqR::merge_kmer_results(resL, resR)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("merge two non empty results that have one k-mer in common", {
  resL <- slam::simple_triplet_matrix(
    i=c(1,3,1,1,2),
    j=c(1,2,3,1,2),
    v=c(1,2,3,4,5),
    nrow=3,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1"))
  )
  
  resR <- slam::simple_triplet_matrix(
    i=c(3,1,1,1,2),
    j=c(1,2,3,1,4),
    v=c(1,2,3,4,5),
    nrow=3,
    dimnames = list(NULL, c("X.X_1", "X.A_1", "X.B.A.A_1.1.0", "A.A_1"))
  )
  
  expected_res <- slam::simple_triplet_matrix(
    i=c(1,3,1,1,2,6,4,4,4,5),
    j=c(1,2,3,1,2,4,5,6,4,1),
    v=c(1,2,3,4,5,1,2,3,4,5),
    nrow=6,
    dimnames = list(NULL, c("A.A_1", "B.A_1", "B.A.A_1.1", "X.X_1", "X.A_1", "X.B.A.A_1.1.0"))
  )
  
  res <- seqR::merge_kmer_results(resL, resR)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})
