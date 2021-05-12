library(testthat)
source("utils.R")

test_that("null alphabet throws an error", {
  expect_error(seqR::count_kmers(sequences=list(c("a", "a", "a")),
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

test_that("alphabet has incompatible element (integer) type with sequences' elements (string list)", {
  expect_error(seqR::count_kmers(sequences=list(c("a", "b")),
                                 alphabet=c(1,2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatible element (numeric) type with sequences' elements (string)", {
  expect_error(seqR::count_kmers(sequences=list(c("aa", "bb")),
                                 alphabet=c(1.2, 2.2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatible element (numeric) type with sequences from vector input (string)", {
  expect_error(seqR::count_kmers(sequences=c("aaaaaaa"),
                                 alphabet=c(1.1, 2.2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatible element (integer) type with sequences from vector input (string)", {
  expect_error(seqR::count_kmers(sequences=c("aaaaaa"),
                                 alphabet=c(1,2),
                                 k=1),
               "alphabet should contain strings")
})

# INVALID K-MER PARAMS ----

test_that("k = 0 generate an error", {
  expect_error(seqR::count_kmers(sequences=c("AAAAA"),
                                 alphabet=c("A"),
                                 k=0),
               "k should be a positive integer")
})

test_that("non integer gaps vector generates an error", {
  expect_error(seqR::count_kmers(sequences=c("AAAAA"),
                                 alphabet=c("A"),
                                 k=1,
                                 kmer_gaps=c("A")),
               "gaps should be an integer vector")
})

test_that("kmer gaps length larger than k-1 generates an error", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                 alphabet=c("A"),
                                 k=1,
                                 kmer_gaps=c(1,2)),
               "the length of kmer_gaps vector should be at most k-1")
})

test_that("unsupported kmer dictionary type raises an error", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=2,
                                kmer_dictionary_name = "unknown"),
               "unsupported k-mer dictionary name: unknown")
})

# BATCH SIZE ----

test_that("provided batch size param is a negative integer", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = -2),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is a positive non integer", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = 2.2),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is zero", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = 0),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is a string", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = "aaaa"),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is an integer vector", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = c(1,2,3)),
               "batch size field must be a positive integer number")
})

test_that("provided batch size param is NULL", {
  expect_error(seqR::count_kmers(sequences=c("AAAA"),
                                alphabet=c("A"),
                                k=1,
                                batch_size = NULL),
               "batch size field must be a positive integer number")
})

# COUNTING ----

test_that("check the result type", {
  res <- seqR::count_kmers("AAAAAAA")
  expect_s3_class(res, "seqR_simple_triplet_matrix")
})

test_that("test list input sequences for gapped k-mers", {
  sq <- c("AAAA", "AAACA")
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
  sq <-c("AAAAAA", "AAACA")
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
  sq <- c("AAAAA", "AA", "AAAAAAAB", "BBB")
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
  sq <- c("aaaaacbb", "aa")
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
  sq <- c("aaaaacbb", "aa", "bb", "aaa")
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
  sq <- c("aa", "aaaaacbb", "bb", "aaa")
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
  sq <- c("AAAAA", "AA", "AAAAAAAB", "BBB")
  
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

test_that("count positional 1-mers for one-dimensional hash P_a, 1_b (P is hashing prime)", {
  P <- 101
  sq <- c(paste0("B", strrep("C", P - 1), "A"))
  
  expected_res <- matrix(c(1, 1), nrow=1, byrow = TRUE)
  colnames(expected_res) <- c("102_A", "1_B")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = c("A", "B"),
                           k = 1,
                           positional = TRUE,
                           hash_dim = 1)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("count positional gapped 2-mers (gap == 1) for one-dimensional hash 1_B.C_1, 10202_A.A_1", {
  sq <- c(paste0("BCA", strrep("C", 10198), "ACA"))
  
  expected_res <- matrix(c(1, 1), nrow=1, byrow = TRUE)
  colnames(expected_res) <- c("1_B.A_1", "10202_A.A_1")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = c("A", "B"),
                           k = 2,
                           kmer_gaps = c(1),
                           positional = TRUE,
                           hash_dim = 1)

  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("count 3-mers without k-mer names, one sequence", {
  sq <- c("AAAAAAAAA")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = LETTERS,
                           k = 3,
                           with_kmer_names = FALSE)
  
  
  
  expect_true(is.null(res$dimnames[[2]]))
  expect_equal(as.matrix(res), matrix(c(7)))
})

test_that("count 3-mers without k-mer names, multiple sequences", {
  sq <- c("AAAAAAAAA", "ACADSDSA", "AAABBB")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = LETTERS,
                           k = 3,
                           with_kmer_names = FALSE)
  
  expect_true(is.null(res$dimnames[[2]]))
})

test_that("(string vector) count 2-mers with alphabet = all", {
  sq <- c("XXXX", "XAXA", "ABC")
  
  expected_res <- matrix(c(
    3, 0, 0, 0, 0,
    0, 2, 1, 0, 0,
    0, 0, 0, 1, 1
  ), byrow=TRUE, nrow=3)
  colnames(expected_res) <- c("X.X_0", "X.A_0", "A.X_0", "A.B_0", "B.C_0")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = "all",
                           k = 2)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("(string vector) count 2-mers with alphabet = all, batch size = 1", {
  sq <- c("XXXX", "XAXA", "ABC")
  
  expected_res <- matrix(c(
    3, 0, 0, 0, 0,
    0, 2, 1, 0, 0,
    0, 0, 0, 1, 1
  ), byrow=TRUE, nrow=3)
  colnames(expected_res) <- c("X.X_0", "X.A_0", "A.X_0", "A.B_0", "B.C_0")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = "all",
                           batch_size = 1,
                           k = 2)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("(string list) count 2-mers with alphabet = all", {
  sq <- list(c("X", "X", "X", "X"),
             c("X", "A", "X", "A"),
             c("A", "B", "C"))
  
  expected_res <- matrix(c(
    3, 0, 0, 0, 0,
    0, 2, 1, 0, 0,
    0, 0, 0, 1, 1
  ), byrow=TRUE, nrow=3)
  colnames(expected_res) <- c("X.X_0", "X.A_0", "A.X_0", "A.B_0", "B.C_0")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = "all",
                           k = 2)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})

test_that("(string list) count 2-mers with alphabet = all, batch_size = 1", {
  sq <- list(c("X", "X", "X", "X"),
             c("X", "A", "X", "A"),
             c("A", "B", "C"))
  
  expected_res <- matrix(c(
    3, 0, 0, 0, 0,
    0, 2, 1, 0, 0,
    0, 0, 0, 1, 1
  ), byrow=TRUE, nrow=3)
  colnames(expected_res) <- c("X.X_0", "X.A_0", "A.X_0", "A.B_0", "B.C_0")
  
  res <- seqR::count_kmers(sequences = sq,
                           alphabet = "all",
                           batch_size = 1,
                           k = 2)
  
  expect_matrices_equal(as.matrix(res), expected_res)
})
