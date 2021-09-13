library(testthat)
source("utils.R")

invoke_test <- function(expected_res,
                        ...) {
  triplet_matrix_res <- seqR::count_multimers(hash_dim=2,
                                              verbose=FALSE,
                                              ...)
  expect_matrices_equal(expected_res, as.matrix(triplet_matrix_res))
}

test_that("(string vector) count 2-mers and 3-mers for list input (AC){100}", {
  expected_res <- matrix(c(
    100, 99, 99, 99
  ), nrow=1, byrow=TRUE)
  colnames(expected_res) <- c("A.C_0", "C.A_0", "A.C.A_0.0", "C.A.C_0.0")
  
  sq <- c(strrep("AC", 100))
  
  invoke_test(expected_res=expected_res,
              sequences=sq,
              k_vector=c(2, 3),
              kmer_alphabet=c("A", "C"),
              positional_vector=rep(FALSE, 2),
              kmer_gaps_list=list(c(), c()),
              with_kmer_counts = TRUE
  )
})

test_that("(string vector) count 2-mers and 3-mers for 2 sequences: (AC){100}, (AD){10}", {
  expected_res <- matrix(c(
    100, 99, 99, 99, 0, 0, 0, 0,
    0,    0,  0,  0, 10, 9, 9, 9
  ), nrow=2, byrow=TRUE)
  colnames(expected_res) <- c("A.C_0", "C.A_0", "A.C.A_0.0", "C.A.C_0.0", "A.D_0", "D.A_0", "A.D.A_0.0", "D.A.D_0.0")
  
  sq <- c(strrep("AC", 100), strrep("AD", 10))
  
  invoke_test(expected_res=expected_res,
              sequences=sq,
              k_vector=c(2, 3),
              kmer_alphabet=c("A", "C", "D"),
              positional_vector=rep(FALSE, 2),
              kmer_gaps_list=list(c(), c()),
              with_kmer_counts = TRUE
  )
})

test_that("(string vector) count 2-mers for 2 sequences: (AC){100}, (AD){10} with defaults", {
  expected_res <- matrix(c(
    100, 99, 99, 99, 0, 0, 0, 0,
    0,    0,  0,  0, 10, 9, 9, 9
  ), nrow=2, byrow=TRUE)
  colnames(expected_res) <- c("A.C_0", "C.A_0", "A.C.A_0.0", "C.A.C_0.0", "A.D_0", "D.A_0", "A.D.A_0.0", "D.A.D_0.0")
  
  sq <- c(strrep("AC", 100), strrep("AD", 10))
  
  res <- seqR::count_multimers(sq, c(2,3), kmer_alphabet=c("A", "D", "C"))
  expect_matrices_equal(expected_res, as.matrix(res))
})

test_that("(string vector) the last sequence do not contain any k-mer", {
  res <- count_multimers(c("aaaaacb", "aa"), c(3, 5), letters)
  expect_equal(nrow(res), 2)
})

test_that("(string vector) expect sparse Matrix as an output", {
  sq <- c(strrep("AC", 100), strrep("AD", 10))
  res <- seqR::count_multimers(sq, c(2,3), kmer_alphabet=c("A", "D", "C"))
  
  expect_is(res, "dgCMatrix")
})

test_that("batch_size = 1 generates message", {
  expect_message(seqR::count_multimers("AAAA", k_vector=c(1, 2)))
})