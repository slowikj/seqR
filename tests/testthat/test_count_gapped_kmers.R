library(testthat)
source("utils.R")

invoke_test <- function(fun, expected_res, ...) {
  res <- seqR::count_kmers(batch_size=200,
                           hash_dim=2,
                           verbose=FALSE,
                           ...)
  expect_matrices_equal(expected_res, as.matrix(res))
}

# STRING LIST TESTS ----

test_that("(string list) test one sequence with kmer_gaps (1) not positional", {
  invoke_test(expected_res=to_matrix(c("a.a_1"=3, "b.b_1"=1, "a.b_1"=1)),
              kmer_alphabet=c("a", "b"),
              sequences=list(c("a", "a", "a", "b", "a", "b", "a")),
              kmer_gaps=c(1),
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string list) test one sequence with kmer_gaps (1) positional", {
  invoke_test(expected_res=to_matrix(c("1_a.a_1"=1, "2_a.b_1"=1, "3_a.a_1"=1, "4_b.b_1"=1, "5_a.a_1"=1)),
              kmer_alphabet=c("a", "b"),
              sequences=list(c("a", "a", "a", "b", "a", "b", "a")),
              kmer_gaps=c(1),
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string list) test one sequence with gapps (1,0) not positional", {
  invoke_test(expected_res=to_matrix(c("b.b.a_1.0"=1, "a.a.b_1.0"=2, "a.b.a_1.0"=1)),
              kmer_alphabet=c("a", "b"),
              sequences=list(c("a", "a", "a", "b", "a", "b", "a")),
              kmer_gaps=c(1, 0),
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string list) test one sequence with gapps (1,0) not positional, alphabet all", {
  invoke_test(expected_res=to_matrix(c("b.b.a_1.0"=1, "a.a.b_1.0"=2, "a.b.a_1.0"=1)),
              kmer_alphabet="all",
              sequences=list(c("a", "a", "a", "b", "a", "b", "a")),
              kmer_gaps=c(1, 0),
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string list) test 2 sequences with kmer_gaps (1,1) positional; some items are not from alphabet", {
  sequences <- list(
    c("a", "b", "c", "c", "as", "b", "c", "a", "a", "b", "a"),
    c("b", "b", "c", "c", "a",  "a", "b", "c", "a", "b", "a"))
  expectedRes <- matrix(c(
    1, 0, 0,
    0, 1, 1
  ), nrow=2, byrow=TRUE)
  colnames(expectedRes) <- c("6_b.a.b_1.1", "7_b.a.a_1.1", "5_a.b.a_1.1")
  
  invoke_test(expected_res=expectedRes,
              kmer_alphabet=c("a", "b"),
              sequences=sequences,
              kmer_gaps=c(1,1),
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string list) the k-mer is longer than a sequence", {
  sequences <- list(
    c("a", "b", "a", "b", "a", "b", "a", "b", "a"),
    c("a", "b", "a", "a", "a", "a", "b", "a", "a"))
  expectedRes <- matrix(nrow=2, ncol=0)
  
  invoke_test(expected_res=expectedRes,
              kmer_alphabet=c("a", "b"),
              sequences=sequences,
              kmer_gaps=rep(1, 10000),
              positional=FALSE,
              with_kmer_counts=TRUE)
})

# STRING VECTOR INPUT TESTS ----

test_that("(string vector) count non positional k-mers (0, 1)", {
  sequences <- c("AAAAAC", "AAA", "AAAC")
  expected_res <- matrix(c(
    2, 1,
    0, 0,
    0, 1), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.1", "A.A.C_0.1")
  invoke_test(expected_res = expected_res,
              kmer_alphabet=c("A", "C"),
              sequences = sequences,
              kmer_gaps = c(0,1),
              positional = FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string vector) count non positional k-mers (0, 1), alphabet all", {
  sequences <- c("AAAAAC", "AAA", "AAAC")
  expected_res <- matrix(c(
    2, 1,
    0, 0,
    0, 1), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.1", "A.A.C_0.1")
  invoke_test(expected_res = expected_res,
              kmer_alphabet="all",
              sequences = sequences,
              kmer_gaps = c(0,1),
              positional = FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string vector) count non positional k-mers (0, 1); some items are not allowed", {
  sequences <- c("AAAACAAAAC", "AACTAAAA", "AACTAAAAC")
  expected_res <- matrix(c(
    3, 0, 0,
    1, 1, 1,
    1, 1, 1
  ), nrow = 3, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.1", "A.A.T_0.1", "T.A.A_0.1")
  invoke_test(expected_res = expected_res,
              kmer_alphabet=c("A", "T"),
              sequences = sequences,
              kmer_gaps = c(0, 1),
              positional = FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string vector) the k-mer is longer than the sequence", {
  sequences <- c("AAAACAAAAC", "AACTAAAA", "AACTAAAAC")
  expected_res <- matrix(nrow=3, ncol=0)
  invoke_test(expected_res = expected_res,
              kmer_alphabet=c("A", "T"),
              sequences = sequences,
              kmer_gaps = rep(1,100),
              positional = FALSE,
              with_kmer_counts=TRUE)
})
