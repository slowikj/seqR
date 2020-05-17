library(testthat)
source("utils.R")

invoke_test <- function(test_fun, expectedRes, ...) {
  result_list <- test_fun(...)
  res <- convert_seqR_list_to_matrix(result_list)
  expect_matrices_equal(res, expectedRes)
}

invoke_test_string <- function(expectedRes, ...) {
  invoke_test(test_fun=seqR::find_kmers_string,
              expectedRes=expectedRes,
              ...)
}

invoke_test_integer <- function(expectedRes, ...) {
  invoke_test(test_fun=seqR::find_kmers_integer,
              expectedRes=expectedRes,
              ...)
}

test_that("count non-positional 2-mers for 2 sequences", {
  sequenceMatrix <- matrix(
    c("a", "a", "b", "a", "a", "b",
      "b", "a", "b", "b", "b", "a"),
    byrow = TRUE,
    nrow = 2
  )
  expectedRes <- matrix(
    c(2, 2, 1, 0,
      0, 1, 2, 2),
    byrow=TRUE,
    nrow=2
  )
  colnames(expectedRes) <- c("a.a_0", "a.b_0", "b.a_0", "b.b_0")
  invoke_test_string(expectedRes=expectedRes,
                     alphabet=c("a", "b"),
                     sequenceMatrix=sequenceMatrix,
                     k=2,
                     positionalKMers=FALSE,
                     withKMerCounts=TRUE,
                     kmerDictionaryName = "unordered_map",
                     batchSize = 100
  )
})

test_that("count positional 2-mers for 2 sequences", {
  sequenceMatrix <- matrix(
    c("a", "a", "b", "a", "a", "b",
      "b", "a", "b", "b", "b", "a"),
    byrow = TRUE,
    nrow = 2
  )
  expectedRes <- matrix(c(
   1, 0, 1, 1, 0, 1, 0, 1, 0,
   0, 1, 1, 0, 1, 0, 1, 0, 1
  ), byrow=TRUE, nrow=2)
  colnames(expectedRes) <- c("1_a.a_0", "1_b.a_0", "2_a.b_0", "3_b.a_0", "3_b.b_0", "4_a.a_0", "4_b.b_0", "5_a.b_0", "5_b.a_0")
  invoke_test_string(expectedRes=expectedRes,
                     alphabet=c("a", "b"),
                     sequenceMatrix=sequenceMatrix,
                     k=2,
                     positionalKMers=TRUE,
                     withKMerCounts=TRUE,
                     kmerDictionaryName = "unordered_map",
                     batchSize = 100)
})

test_that("count non positional 5-mers for 10 sequences (10^6 each)", {
  nrow <- 10
  ncol <- 1000000
  sequenceMatrix <- matrix(rep(rep(5, ncol), nrow), byrow=TRUE, nrow=nrow)
  
  expectedRes <- matrix(rep(999996, nrow),
                        nrow=nrow)
  colnames(expectedRes) <- c("5.5.5.5.5_0.0.0.0")
  invoke_test_integer(expectedRes=expectedRes,
                      alphabet=c(5),
                      sequenceMatrix=sequenceMatrix,
                      k=5,
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 100)
})

test_that("find non positional 5-mers for 100 sequences (10^6 each) (without k-mer counts)", {
  nrow <- 100
  ncol <- 1000000
  sequenceMatrix <- matrix(rep(rep(5, ncol), nrow), byrow=TRUE, nrow=nrow)

  expectedRes <- matrix(rep(1, nrow),
                        nrow=nrow)
  colnames(expectedRes) <- c("5.5.5.5.5_0.0.0.0")
  invoke_test_integer(expectedRes=expectedRes,
                      alphabet=c(5),
                      sequenceMatrix=sequenceMatrix,
                      k=5,
                      positionalKMers=FALSE,
                      withKMerCounts=FALSE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 100)
})
