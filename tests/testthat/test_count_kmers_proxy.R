library(testthat)

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

test_that("alphabet has incompatible element (numeric) type with sequences from tidysq (string)", {
  expect_error(seqR::count_kmers(sequences=tidysq::as.sq("aaaaaaa"),
                                 alphabet=c(1.1, 2.2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("alphabet has incompatible element (integer) type with sequences from tidysq (string)", {
  expect_error(seqR::count_kmers(sequences=tidysq::as.sq("aaaaaa"),
                                 alphabet=c(1,2),
                                 k=1),
               "alphabet should contain strings")
})

test_that("sequences of unsupported type generate an error", {
  expect_error(seqR::count_kmers(sequences=list(1,2,3),
                                 alphabet=c(1,2,3),
                                 k=1),
               "sequences param has unsupported type")
})

test_that("test tidysq for gapped k-mers", {
  sq <- tidysq::as.sq(c("AAAA", "AAACA"))
  expected_res <- matrix(c(
    2, 0,
    2, 1), byrow=TRUE, nrow=2)
  colnames(expected_res) <- c("A.A", "A.C")
  
  res <- seqR::count_kmers(sequences=sq,
                           k=2,
                           alphabet=c("A", "C"),
                           positional=FALSE,
                           kmer_gaps=c(1))
  expect_equal(expected_res, res)
})
