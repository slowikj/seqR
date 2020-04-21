library(testthat)
source("utils.R")

invoke_test <- function(fun, expected_res, alphabet, sequence, k, positionalKMers) {
  sequenceMatrix <- to_matrix(sequence)
  expected_res <- to_matrix(expected_res)
  
  res <- fun(alphabet=alphabet,
             sequenceMatrix = sequenceMatrix,
             k=k,
             positionalKMers = positionalKMers)
  
  expect_matrices_equal(res, expected_res)
}

invoke_test_string <- function(...) {
  invoke_test(seqR::count_kmers_string, ...)
}

invoke_test_integer <- function(...) {
  invoke_test(seqR::count_kmers_integer, ...)
}

invoke_test_numeric <- function(...) {
  invoke_test(seqR::count_kmers_numeric, ...)
}

# STRING MATRIX TESTS ----

test_that("(string) count non positional 2-mers", {
  invoke_test_string(expected_res = c("a.a_0"=1, "a.b_0"=2, "b.a_0"=2),
                     alphabet=c("a", "b"),
                     sequence=c("a", "b", "a", "b", "a", "a"),
                     k=2,
                     positionalKMers=FALSE)
})

test_that("(string) count positional 2-mers", {
  invoke_test_string(expected_res = c("1_a.b_0"=1, "2_b.a_0"=1, "3_a.b_0"=1, "4_b.a_0"=1, "5_a.a_0"=1),
                     alphabet=c("a", "b"),
                     sequence=c("a", "b", "a", "b", "a", "a"),
                     k=2,
                     positionalKMers=TRUE)
})

test_that("(string) count non positional 1-mers", {
  invoke_test_string(expected_res=c("a"=3, "b"=2),
                     alphabet=c("a", "b"),
                     sequence=c("a", "a", "b", "a", "b"),
                     k=1,
                     positionalKMers=FALSE)
})

test_that("(string) count positional 1-mers", {
  invoke_test_string(expected_res=c("1_a"=1, "2_a"=1, "3_b"=1, "4_a"=1, "5_b"=1),
                     alphabet=c("a", "b"),
                     sequence=c("a", "a", "b", "a", "b"),
                     k=1,
                     positionalKMers=TRUE)
})

test_that("(string) count non positional 1-mers if some sequence items are not allowed", {
  invoke_test_string(expected_res=c("a"=3),
                     alphabet=c("a"),
                     sequence=c("a", "a", "b", "a", "b"),
                     k=1,
                     positionalKMers=FALSE)
})

test_that("(string) count positional 1-mers if some sequence items are not allowed", {
  invoke_test_string(expected_res=c("1_a"=1, "2_a"=1, "4_a"=1),
                     alphabet=c("a"),
                     sequence=c("a", "a", "b", "a", "b"),
                     k=1,
                     positionalKMers=TRUE)
})

test_that("(string) count non positional 3-mers", {
  invoke_test_string(expected_res=c("a.a.b_0.0"=2, "a.b.c_0.0"=2, "b.c.a_0.0"=2, "c.a.a_0.0"=1),
                     alphabet=c("a", "b", "c"),
                     sequence=c("a", "a", "b", "c", "a", "a", "b", "c", "a"),
                     k=3,
                     positionalKMers=FALSE)
})

test_that("(string) count positional 3-mers", {
  invoke_test_string(expected_res=c("1_a.a.a_0.0"=1, "2_a.a.b_0.0"=1, "3_a.b.a_0.0"=1, "4_b.a.a_0.0"=1, "5_a.a.a_0.0"=1),
                     alphabet=c("a", "b"),
                     sequence=c("a", "a", "a", "b", "a", "a", "a"),
                     k=3,
                     positionalKMers=TRUE)
})

test_that("(string) count non positional 1-mers", {
  invoke_test_string(expected_res=c("a"=5, "b"=3, "c"=2),
                     alphabet=c("a", "b", "c"),
                     sequence=c("a", "a", "b", "b", "a", "a", "b", "c", "c", "a"),
                     k=1,
                     positionalKMers=FALSE)
})

test_that("(string) count 2-mers: long (ab){1000} sequence", {
  invoke_test_string(expected_res=c("a.b_0"=1000, "b.a_0"=999),
                     alphabet=c("a", "b"),
                     sequence=rep(c("a", "b"), 1000),
                     k=2,
                     positionalKMers=FALSE)
})

# INTEGER MATRIX TESTS ----

test_that("(integer) count non positional 2-mers with not allowed item", {
  invoke_test_integer(expected_res=c("1.0_0"=2, "1.1_0"=3, "0.1_0"=2),
                      alphabet=c(0, 1),
                      sequence=c(1,0,1,1,0,1,1,2,1,1),
                      k=2,
                      positionalKMers=FALSE)
})

# NUMERIC MATRIX TESTS ----

test_that("(numeric) count non positional 2-mers with not allowed item", {
  invoke_test_numeric(expected_res=c("0.000.1.000_0"=2, "1.000.1.000_0"=3, "1.000.0.000_0"=2),
                     alphabet=c(0, 1),
                     sequence=as.numeric(c(1,0,1,1,0,1,1,2,1,1)),
                     k=2,
                     positionalKMers=FALSE)
})
