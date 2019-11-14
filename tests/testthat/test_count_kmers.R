library(testthat)

test_that("count positional unigrams", {
  res <- count_kmers(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c(),
    c("a", "b", "c"),
    pos=TRUE)
  
  expected <- c("0_c"=1, "2_a"=1, "1_b"=2, "0_a"=1, "2_c"=1)
  expect_mapequal(res, expected)
})

test_that("count non positional unigrams", {
  res <- count_kmers(
    matrix(data=c("a", "c", "b", "b", "a", "c"), nrow=2),
    c(),
    c("a", "b", "c"),
    pos=FALSE)
  
  expected <- c("a"=2, "b"=2, "c"=2)
  expect_mapequal(res, expected)
})

test_that("count non positional 2grams", {
  res <- count_kmers(
    matrix(data=c("a", "b", "c", "a", "b", "c"), nrow=2),
    c(0),
    c("a", "b", "c"),
    FALSE)
  expected <- c("c.c"=0, "b.b"=0, "b.a"=1, "a.a"=0, "c.a"=0, "b.c"=0, "a.c"=2, "a.b"=0, "c.b"=1)
  expect_mapequal(res, expected)
})

test_that("count positional 2grams", {
  res <- count_kmers(
    matrix(data=c("a", "b", "c", "a", "b", "c"), nrow=2),
    c(0),
    c("a", "b", "c"),
    TRUE)
  expected <- c("0_b.a"=1, "0_a.c"=1, "1_a.c"=1, "1_c.b"=1)
  expect_mapequal(res, expected)
})

test_that("count positional 2grams with items longer than 1", {
  res <- count_kmers(
    matrix(data=c("a", "bb", "ss", "aa", "a", "b", "cc", "aa", "a"), nrow=3),
    c(0),
    c("a", "bb", "ss", "aa", "b", "cc"),
    TRUE)
  expected <- c("1_b.a"=1, "1_aa.cc"=1,  "0_bb.a"=1, "0_a.aa"=1,  "0_ss.b"=1, "1_a.aa"=1)
  expect_mapequal(res, expected)
})

test_that("count non positional 2grams with items longer than 1", {
  res <- count_kmers(
    matrix(data=c("a", "bb", "ss", "aa", "a", "b", "cc", "aa", "a"), nrow=3),
    c(0),
    c("a", "bb", "ss", "aa", "b", "cc"),
    FALSE)
  expected <- c("ss.ss"=0, "a.a"=0, "cc.bb"=0,  "b.aa"=0, "cc.aa"=0,  "b.cc"=0, "cc.b"=0,
                "aa.bb"=0,  "b.ss"=0, "ss.a"=0, "ss.bb"=0, "a.aa"=2, "ss.aa"=0, "ss.b"=1,
                "a.bb"=0, "aa.ss"=0,  "cc.a"=0,  "aa.a"=0,  "bb.b"=0, "cc.cc"=0,
                "bb.a"=1, "aa.aa"=0,  "a.cc"=0,  "aa.b"=0,  "a.ss"=0, "ss.cc"=0,  "b.bb"=0,
                "b.a"=1, "bb.aa"=0, "bb.bb"=0, "bb.ss"=0,  "aa.cc"=1, "a.b"=0, "bb.cc"=0, "b.b"=0, "cc.ss"=0)
  expect_mapequal(res, expected)
})
