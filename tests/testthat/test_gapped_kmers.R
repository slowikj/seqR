library(testthat)
source("utils.R")

invoke_test <- function(fun, expected_res, ...) {
  result_list <- fun(...)
  res <- convert_seqR_list_to_matrix(result_list)
  expect_matrices_equal(expected_res, res)
}

invoke_test_strings <- function(...) {
  invoke_test(seqR::find_gapped_kmers_string, ...)
}

invoke_test_integer <- function(...) {
  invoke_test(seqR::find_gapped_kmers_integer, ...)
}

invoke_test_numeric <- function(...) {
  invoke_test(seqR::find_gapped_kmers_numeric, ...)
}

invoke_test_list <- function(...) {
  invoke_test(seqR::find_gapped_kmers_list, ...)
}

# STRING MATRIX TESTS ----

test_that("(string)test one sequence with gaps (1) not positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test_strings(expected_res=to_matrix(c("a.a_1"=3, "b.b_1"=1, "a.b_1"=1)),
                      alphabet=c("a", "b"),
                      sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
                      gaps=c(1),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

test_that("(string)test one sequence with gaps (1) positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test_strings(expected_res=to_matrix(c("1_a.a_1"=1, "2_a.b_1"=1, "3_a.a_1"=1, "4_b.b_1"=1, "5_a.a_1"=1)),
                      alphabet=c("a", "b"),
                      sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
                      gaps=c(1),
                      positionalKMers=TRUE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

test_that("(string)test one sequence with gapps (1,0) not positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test_strings(expected_res=to_matrix(c("b.b.a_1.0"=1, "a.a.b_1.0"=2, "a.b.a_1.0"=1)),
                      alphabet=c("a", "b"),
                      sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
                      gaps=c(1, 0),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

test_that("(string)test 2 sequences with gaps (1,1) positional; some items are not from alphabet", {
  seqMatrix <- matrix(c(
    "a", "b", "c", "c", "as", "b", "c", "a", "a", "b", "a",
    "b", "b", "c", "c", "a",  "a", "b", "c", "a", "b", "a"
  ), nrow=2, byrow=TRUE)
  expectedRes <- matrix(c(
    1, 0, 0,
    0, 1, 1
  ), nrow=2, byrow=TRUE)
  colnames(expectedRes) <- c("6_b.a.b_1.1", "7_b.a.a_1.1", "5_a.b.a_1.1")
  
  invoke_test_strings(expected_res=expectedRes,
                      alphabet=c("a", "b"),
                      sequenceMatrix=seqMatrix,
                      gaps=c(1,1),
                      positionalKMers=TRUE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

test_that("(string) the k-mer is longer than a sequence", {
  seqMatrix <- matrix(c(
    "a", "b", "a", "b", "a", "b", "a", "b", "a",
    "a", "b", "a", "a", "a", "a", "b", "a", "a"
  ), nrow=2, byrow=TRUE)
  expectedRes <- matrix(nrow=2, ncol=0)
  
  invoke_test_strings(expected_res=expectedRes,
                      alphabet=c("a", "b"),
                      sequenceMatrix=seqMatrix,
                      gaps=rep(1, 10000),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

# INTEGER MATRIX TESTS ----

test_that("(integer) test 2 sequences with gaps (1,1) non positional; some items are not from alphabet", {
  seqMatrix <- matrix(c(
    1, 2, 0, 1, 1, 1, 2, 1, 0, 1,
    1, 1, 1, 0, 0, 0, 0, 2, 2, 1
  ), nrow=2, byrow=TRUE)
  expected_res <- matrix(c(
    1, 2, 0, 0,
    0, 0, 1, 2 
  ), nrow=2, byrow=TRUE)
  colnames(expected_res) <- c("1.0.1_1.1", "1.1.1_1.1", "1.1.0_1.1", "1.0.0_1.1")
  invoke_test_integer(expected_res=expected_res,
                      alphabet=c(0,1),
                      sequenceMatrix=seqMatrix,
                      gaps=c(1,1),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

test_that("(integer) the k-mer is longer than the sequence", {
  seqMatrix <- matrix(c(
    1,2,3,2,2,3,2,3,23,2,3,2,3,2,3,2,
    1,2,2,2,3,2,3,2,2,22,2,2,2,2,2,2
  ), nrow=2, byrow=TRUE)
  expectedRes <- matrix(nrow=2, ncol=0)
  invoke_test_integer(expected_res=expectedRes,
                      alphabet=1:10,
                      sequenceMatrix=seqMatrix,
                      gaps=rep(1, 100),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

# NUMERIC MATRIX TESTS ----

test_that("(numeric) test 2 sequences with gaps (1,1) non positional; some items are not from alphabet", {
  seqMatrix <- matrix(as.numeric(c(
    1, 2, 0, 1, 1, 1, 2, 1, 0, 1,
    1, 1, 1, 0, 0, 0, 0, 2, 2, 1
  )), nrow=2, byrow=TRUE)
  expected_res <- matrix(c(
    1, 2, 0, 0,
    0, 0, 1, 2 
  ), nrow=2, byrow=TRUE)
  colnames(expected_res) <- c("1.000.0.000.1.000_1.1", "1.000.1.000.1.000_1.1", "1.000.1.000.0.000_1.1", "1.000.0.000.0.000_1.1")
  invoke_test_numeric(expected_res=expected_res,
                      alphabet=c(0,1),
                      sequenceMatrix=seqMatrix,
                      gaps=c(1,1),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

test_that("(numeric) the k-mer is longer than the sequence", {
  seqMatrix <- matrix(as.numeric(c(
    1,2,3,2,2,3,2,3,23,2,3,2,3,2,3,2,
    1,2,2,2,3,2,3,2,2,22,2,2,2,2,2,2
  )), nrow=2, byrow=TRUE)
  expectedRes <- matrix(nrow=2, ncol=0)
  invoke_test_numeric(expected_res=expectedRes,
                      alphabet=as.numeric(1:10),
                      sequenceMatrix=seqMatrix,
                      gaps=rep(1, 100),
                      positionalKMers=FALSE,
                      withKMerCounts=TRUE,
                      kmerDictionaryName = "unordered_map",
                      batchSize = 200)
})

# LIST INPUT TESTS ----

test_that("(list input) count non positional k-mers (0, 1)", {
  sq <- list("AAAAAC", "AAA", "AAAC")
  expected_res <- matrix(c(
    2, 1,
    0, 0,
    0, 1), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.1", "A.A.C_0.1")
  invoke_test_list(expected_res = expected_res,
                      alphabet=c("A", "C"),
                      sq = sq,
                      gaps = c(0,1),
                      positionalKMers = FALSE,
                     withKMerCounts=TRUE,
                     kmerDictionaryName = "unordered_map",
                     batchSize = 200)
})

test_that("(list input) count non positional k-mers (0, 1); some items are not allowed", {
  sq <- list("AAAACAAAAC", "AACTAAAA", "AACTAAAAC")
  expected_res <- matrix(c(
    3, 0, 0,
    1, 1, 1,
    1, 1, 1
  ), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.1", "A.A.T_0.1", "T.A.A_0.1")
  invoke_test_list(expected_res = expected_res,
                     alphabet=c("A", "T"),
                     sq = sq,
                     gaps = c(0, 1),
                     positionalKMers = FALSE,
                     withKMerCounts=TRUE,
                     kmerDictionaryName = "unordered_map",
                     batchSize = 200)
})

test_that("(list input) the k-mer is longer than the sequence", {
  sq <- list("AAAACAAAAC", "AACTAAAA", "AACTAAAAC")
  expected_res <- matrix(nrow=3, ncol=0)
  invoke_test_list(expected_res = expected_res,
                     alphabet=c("A", "T"),
                     sq = sq,
                     gaps = rep(1,100),
                     positionalKMers = FALSE,
                     withKMerCounts=TRUE,
                     kmerDictionaryName = "unordered_map",
                     batchSize = 200)
})

