library(testthat)
source("utils.R")

invoke_test <- function(fun, expected_res, seq, ...) {
  expected_res <- to_matrix(expected_res)
  
  res <- seqR::count_kmers(sequences=seq,
                           batch_size=200,
                           hash_dim=2,
                           verbose=FALSE,
                           ...)
  
  expect_matrices_equal(as.matrix(res), expected_res)
}

# STRING LIST TESTS ----

test_that("(string) count non positional 2-mers", {
  invoke_test(expected_res = c("a.a_0"=1, "a.b_0"=2, "b.a_0"=2),
              kmer_alphabet=c("a", "b"),
              seq=list(c("a", "b", "a", "b", "a", "a")),
              k=2,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count positional 2-mers", {
  invoke_test(expected_res = c("1_a.b_0"=1, "2_b.a_0"=1, "3_a.b_0"=1, "4_b.a_0"=1, "5_a.a_0"=1),
              kmer_alphabet=c("a", "b"),
              seq=list(c("a", "b", "a", "b", "a", "a")),
              k=2,
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 1-mers", {
  invoke_test(expected_res=c("a"=3, "b"=2),
              kmer_alphabet=c("a", "b"),
              seq=list(c("a", "a", "b", "a", "b")),
              k=1,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count positional 1-mers", {
  invoke_test(expected_res=c("1_a"=1, "2_a"=1, "3_b"=1, "4_a"=1, "5_b"=1),
              kmer_alphabet=c("a", "b"),
              seq=list(c("a", "a", "b", "a", "b")),
              k=1,
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 1-mers if some sequence items are not allowed", {
  invoke_test(expected_res=c("a"=3),
              kmer_alphabet=c("a"),
              seq=list(c("a", "a", "b", "a", "b")),
              k=1,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count positional 1-mers if some sequence items are not allowed", {
  invoke_test(expected_res=c("1_a"=1, "2_a"=1, "4_a"=1),
              kmer_alphabet=c("a"),
              seq=list(c("a", "a", "b", "a", "b")),
              k=1,
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 3-mers", {
  invoke_test(expected_res=c("a.a.b_0.0"=2, "a.b.c_0.0"=2, "b.c.a_0.0"=2, "c.a.a_0.0"=1),
              kmer_alphabet=c("a", "b", "c"),
              seq=list(c("a", "a", "b", "c", "a", "a", "b", "c", "a")),
              k=3,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count positional 3-mers", {
  invoke_test(expected_res=c("1_a.a.a_0.0"=1, "2_a.a.b_0.0"=1, "3_a.b.a_0.0"=1, "4_b.a.a_0.0"=1, "5_a.a.a_0.0"=1),
              kmer_alphabet=c("a", "b"),
              seq=list(c("a", "a", "a", "b", "a", "a", "a")),
              k=3,
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 1-mers", {
  invoke_test(expected_res=c("a"=5, "b"=3, "c"=2),
              kmer_alphabet=c("a", "b", "c"),
              seq=list(c("a", "a", "b", "b", "a", "a", "b", "c", "c", "a")),
              k=1,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count 2-mers: long (ab){1000} sequence", {
  invoke_test(expected_res=c("a.b_0"=1000, "b.a_0"=999),
              kmer_alphabet=c("a", "b"),
              seq=list(rep(c("a", "b"), 1000)),
              k=2,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 2-mers which contains only 'a' character", {
  invoke_test(expected_res=c("a.a_0"=4),
              kmer_alphabet=c("a"),
              seq=list(c("a", "b", "c", "a", "a", "b", "a", "a", "a", "a", "b")),
              k=2,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count positional 2-mers which contains only 'a' character", {
  invoke_test(expected_res=c("4_a.a_0"=1, "7_a.a_0"=1, "8_a.a_0"=1, "9_a.a_0"=1),
              kmer_alphabet=c("a"),
              seq=list(c("a", "b", "c", "a", "a", "b", "a", "a", "a", "a", "b")),
              k=2,
              positional=TRUE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 2-mers which contains only 'a' or 'b' characters", {
  invoke_test(expected_res=c("a.a_0"=1, "a.b_0"=2, "b.a_0"=1),
              kmer_alphabet=c("a", "b"),
              seq=list(c("x", "x", "a", "x", "b", "x", "x", "a", "a", "b", "a", "aa", "a", "b", "x", "a")),
              k=2,
              positional=FALSE,
              with_kmer_counts=TRUE)
})

test_that("(string) count non positional 2-mers which contains only 'a', 'b', 'x' characters", {
  invoke_test(expected_res=c("a.a_0"=1, "a.b_0"=2, "b.a_0"=1, "x.x_0"=2, "x.a_0"=3, "a.x_0"=1, "b.x_0"=2, "x.b_0"=1),
              kmer_alphabet=c("a", "b", "x"),
              seq=list(c("x", "x", "a", "x", "b", "x", "x", "a", "a", "b", "a", "aa", "a", "b", "x", "a")),
              k=2,
              positional=FALSE,
              with_kmer_counts=TRUE)
})
