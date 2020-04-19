library(testthat)
source("utils.R")

invoke_test <- function(fun, expected_res, ...) {
  res <- fun(...)
  expect_matrices_equal(expected_res, res)
}

invoke_test_strings <- function(...) {
  invoke_test(seqR::count_gapped_kmers_string, ...)
}

invoke_test_integer <- function(...) {
  invoke_test(seqR::count_gapped_kmers_integer, ...)
}

invoke_test_numeric <- function(...) {
  invoke_test(seqR::count_gapped_kmers_numeric, ...)
}

invoke_test_tidysq <- function(...) {
  invoke_test(seqR::count_gapped_kmers_tidysq, ...)
}

# STRING

test_that("(string)test one sequence with gaps (1) not positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test_strings(expected_res=to_matrix(c("a.a_1"=3, "b.b_1"=1, "a.b_1"=1)),
                      alphabet=c("a", "b"),
                      sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
                      gaps=c(1),
                      positionalKMers=FALSE)
})

test_that("(string)test one sequence with gaps (1) positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test_strings(expected_res=to_matrix(c("1_a.a_1"=1, "2_a.b_1"=1, "3_a.a_1"=1, "4_b.b_1"=1, "5_a.a_1"=1)),
                      alphabet=c("a", "b"),
                      sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
                      gaps=c(1),
                      positionalKMers=TRUE)
})

test_that("(string)test one sequence with gapps (1,0) not positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test_strings(expected_res=to_matrix(c("b.b.a_1.0"=1, "a.a.b_1.0"=2, "a.b.a_1.0"=1)),
                      alphabet=c("a", "b"),
                      sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
                      gaps=c(1, 0),
                      positionalKMers=FALSE)
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
                      positionalKMers=TRUE)
})

# INTEGER

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
                      positionalKMers=FALSE)
})

# NUMERIC

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
                      positionalKMers=FALSE)
})

# TIDYSQ

test_that("(tidysq) count non positional k-mers (0, 1)", {
  sq <- tidysq::construct_sq(c("AAAAAC", "AAA", "AAAC"), type="ami")
  expected_res <- matrix(c(
    2, 1,
    0, 0,
    0, 1), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.1", "A.A.C_0.1")
  invoke_test_tidysq(expected_res = expected_res,
                      alphabet=c("A", "C"),
                      sq = sq,
                      gaps = c(0,1),
                      positionalKMers = FALSE)
})

