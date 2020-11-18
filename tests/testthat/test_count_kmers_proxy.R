library(testthat)
source("utils.R")

test_that("null alphabet throws an error", {
  expect_error(seqR::count_kmers(sequences=c("a", "a", "a"),
                                alphabet=c()),
               "alphabet param is empty")
})

test_that("null sequences throws an error", {
  expect_error(seqR::count_kmers(sequences=c(),
                                 alphabet=c("a"),
                                 k=1),
               "sequences param is empty")
})

# ALPHABET ----

test_that("alphabet has incompatible element (integer) type with sequences' elements (string)", {
  expect_error(seqR::count_kmers(sequences=c("a", "b"),
                                 alphabet=c(1,2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatible element (string) type with sequences' elements (integer)", {
  expect_error(seqR::count_kmers(sequences=c(1,2),
                                 alphabet=c("a", "b"),
                                 k=1),
               "alphabet should contain integers")
})

test_that("alphabet has incompatible element (numeric) type with sequences' elements (integer)", {
  expect_error(seqR::count_kmers(sequences=c(1, 2),
                                 alphabet=c(1.5, 2.2),
                                 k=1),
               "alphabet should contain integers")
})

test_that("alphabet has incompatible element (integer) type with sequences' elements (numeric)", {
  expect_error(seqR::count_kmers(sequences=c(1.1, 2.2),
                                 alphabet=c(1, 2),
                                 k=1),
               "alphabet should contain numerics")
})

test_that("alphabet has incompatible element (numeric) type with sequences' elements (string)", {
  expect_error(seqR::count_kmers(sequences=c("aa", "bb"),
                                 alphabet=c(1.2, 2.2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatile element (string) type with sequences' elements  (numeric)", {
  expect_error(seqR::count_kmers(sequences=c(1.2, 1.1),
                                 alphabet=c("aa", "bb"),
                                 k=1),
               "alphabet should contain numerics")
})

test_that("alphabet has incompatible element (numeric) type with sequences from list input (string)", {
  expect_error(seqR::count_kmers(sequences=list("aaaaaaa"),
                                 alphabet=c(1.1, 2.2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatible element (integer) type with sequences from list input (string)", {
  expect_error(seqR::count_kmers(sequences=list("aaaaaa"),
                                 alphabet=c(1,2),
                                 k=1),
               "alphabet should contain strings")
})

# SEQUENCE TYPE ----

# INVALID K-MER PARAMS ----

test_that("k = 0 generate an error", {
  expect_error(seqR::count_kmers(sequences=list("AAAAA"),
                                 alphabet=c("A"),
                                 k=0),
               "k should be a positive integer")
})

test_that("non integer gaps vector generates an error", {
  expect_error(seqR::count_kmers(sequences=list("AAAAA"),
                                 alphabet=c("A"),
                                 k=1,
                                 kmer_gaps=c("A")),
               "gaps should be an integer vector")
})

test_that("kmer gaps length larger than k-1 generates an error", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                 alphabet=c("A"),
                                 k=1,
                                 kmer_gaps=c(1,2)),
               "the length of kmer_gaps vector should be at most k-1")
})

test_that("unsupported kmer dictionary type raises an error", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=2,
                                kmer_dictionary_name = "unknown"),
               "unsupported k-mer dictionary name: unknown")
})

# BATCH SIZE ----

test_that("provided batch size param is a negative integer", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = -2),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is a positive non integer", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = 2.2),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is zero", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = 0),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is a string", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = "aaaa"),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is an integer vector", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = c(1,2,3)),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is NULL", {
  expect_error(seqR::count_kmers(sequences=list("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = NULL),
               "batch size field must be a positive integer number")
})

# COUNTING ----

test_that("test list input sequences for gapped k-mers", {
  sq <- list("AAAA", "AAACA")
  expected_res <- matrix(c(
    2, 0,
    2, 1), byrow=TRUE, nrow=2)
  colnames(expected_res) <- c("A.A_1", "A.C_1")
  
  res <- seqR::count_kmers(sequences=sq,
                          k=2,
                          alphabet=c("A", "C"),
                          positional=FALSE,
                          kmer_gaps=c(1),
                          kmer_dictionary_name = "unordered_map")
  
  expect_equal(expected_res, as.matrix(res))
})

test_that("test the case when there is 0 found k-mers", {
  sq <-list("AAAAAA", "AAACA")
  expected_res <- matrix(nrow=2, ncol=0)
  res <- seqR::count_kmers(sequences=sq,
                          k=100,
                          alphabet=c("A"),
                          positional=FALSE,
                          kmer_gaps=c(1),
                          kmer_dictionary_name="unordered_map")
  
  expect_equal(expected_res, as.matrix(res))
})

# DIFFERENT BATCHES SIZES ----

run_batch_test <- function(batch_size) {
  sq <- list("AAAAA", "AA", "AAAAAAAB", "BBB")
  expected_res <- matrix(c(
    3, 0, 0,
    0, 0, 0,
    5, 1, 0,
    0, 0, 1
  ), nrow=4, byrow=TRUE)
  colnames(expected_res) <- c("A.A.A_0.0", "A.A.B_0.0", "B.B.B_0.0")
  res <- seqR::count_kmers(sequences = sq,
                          alphabet=c("A", "B"),
                          k=3,
                          positional=FALSE,
                          with_kmer_counts=TRUE,
                          kmer_dictionary_name="linear_list",
                          batch_size = batch_size)
  expect_equal(expected_res, as.matrix(res))
}

test_that("test list input sequences that are processed in ONE batch iteration", {
  run_batch_test(batch_size=3)
})

test_that("test list input sequences that are processed in TWO batch iterations", {
  run_batch_test(batch_size=2)
})

test_that("test list input sequences that are processed in THREE batch iterations", {
  run_batch_test(batch_size=1)
})

# CORNER CASES ----

test_that("the last input sequence does not contain any specified k-mer", {
  sq <- list("aaaaacbb", "aa")
  expected_res <- matrix(c(
    1, 1, 1, 1,
    0, 0, 0, 0
  ), nrow=2, byrow=TRUE)
  colnames(expected_res) <- c("a.a.c.b.b_0.0.0.0", "a.a.a.c.b_0.0.0.0", "a.a.a.a.a_0.0.0.0", "a.a.a.a.c_0.0.0.0")
  
  res <- seqR::count_kmers(sequences=sq,
                           alphabet=letters,
                           k=5,
                           positional=FALSE,
                           with_kmer_counts = FALSE,
                           kmer_dictionary_name = "unordered_map",
                           batch_size = 100)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("more than one last input sequences do not contain any specified k-mer", {
  sq <- list("aaaaacbb", "aa", "bb", "aaa")
  expected_res <- matrix(c(
    1, 1, 1, 1,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0
  ), nrow=4, byrow=TRUE)
  colnames(expected_res) <- c("a.a.c.b.b_0.0.0.0", "a.a.a.c.b_0.0.0.0", "a.a.a.a.a_0.0.0.0", "a.a.a.a.c_0.0.0.0")
  
  res <- seqR::count_kmers(sequences=sq,
                           alphabet=letters,
                           k=5,
                           positional=FALSE,
                           with_kmer_counts = FALSE,
                           kmer_dictionary_name = "unordered_map",
                           batch_size = 100)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("some input sequences do not contain any specified k-mer", {
  sq <- list("aa", "aaaaacbb", "bb", "aaa")
  expected_res <- matrix(c(
    0, 0, 0, 0,
    1, 1, 1, 1,
    0, 0, 0, 0,
    0, 0, 0, 0
  ), nrow=4, byrow=TRUE)
  colnames(expected_res) <- c("a.a.c.b.b_0.0.0.0", "a.a.a.c.b_0.0.0.0", "a.a.a.a.a_0.0.0.0", "a.a.a.a.c_0.0.0.0")
  
  res <- seqR::count_kmers(sequences=sq,
                           alphabet=letters,
                           k=5,
                           positional=FALSE,
                           with_kmer_counts = FALSE,
                           kmer_dictionary_name = "unordered_map",
                           batch_size = 100)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("expect simple_triplet_matrix as an output", {
  sq <- list("AAAAA", "AA", "AAAAAAAB", "BBB")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet=c("A", "B"),
                           k=3,
                           positional=FALSE,
                           with_kmer_counts=TRUE,
                           kmer_dictionary_name="linear_list",
                           batch_size = 100)
  
  expect_is(res, "simple_triplet_matrix")
})

# POSITIONAL K-MERS ----

test_that("count positional 1-mers for one-dimensional hash P_a, 0_b (P is hashing prime)", {
  sq <- list(paste0("B", strrep("C", 100), "A"))
  
  expected_res <- matrix(c(1, 1), nrow=1, byrow = TRUE)
  colnames(expected_res) <- c("102_A", "1_B")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = c("A", "B"),
                           k = 1,
                           positional = TRUE,
                           hash_dim = 1)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})
