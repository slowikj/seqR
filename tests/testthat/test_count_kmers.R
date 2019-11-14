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

test_that("count non positional gapped 3 grams: 1_1_1", {
  res <- count_kmers(
    c("a", "a", "c", "b", "b", "a", "c"),
    c(1,1),
    c("a", "b", "c"),
    FALSE
  )
  
  expected <- c("b.a.b"=0, "c.b.c"=1, "a.a.a"=0, "a.b.a"=1, "a.b.b"=0, "c.a.c"=0,
                "c.b.a"=0, "a.c.c"=0, "b.b.c"=0, "c.b.b"=0, "c.c.a"=0, "c.a.a"=0, "a.a.c"=0,
                "b.b.a"=0, "c.c.b"=0, "b.c.b"=0, "c.c.c"=0, "b.c.a"=0, "b.a.a"=0, "a.b.c"=0,
                "c.a.b"=0, "a.a.b"=0, "a.c.a"=0, "b.b.b"=0, "a.c.b"=1, "b.c.c"=0, "b.a.c"=0)
  expect_mapequal(res, expected)
})

test_that("count positional gapped 3 grams: 1_1_1", {
  res <- count_kmers(
    c("a", "a", "c", "b", "b", "a", "c"),
    c(1,1),
    c("a", "b", "c"),
    TRUE
  )
  
  expected <- c("0_a.c.b"=1, "2_c.b.c"=1, "1_a.b.a"=1)
  expect_mapequal(res, expected)
})

test_that("count positional gapped 4 grams: 1_1__1_1", {
  res <- count_kmers(
    c("a", "b", "c", "aaaa", "123", "5dd", "a", "a", "a"),
    c(1,2,1),
    c("a", "b", "c", "aaaa", "123", "5dd"),
    TRUE
  )
  
  expected <- c("0_a.c.5dd.a"=1, "1_b.aaaa.a.a"=1)
  expect_mapequal(res, expected)
})

test_that("count all non positional 3grams 1_2 in rep('a', 10)", {
  res <- count_kmers(
    rep("a", 10),
    c(1,2),
    c("a"),
    FALSE
  )
  expected <- c("a.a.a"=5)
  expect_mapequal(res, expected)
})

test_that("count all non positional 3 grams with extended alphabet in rep('a', 10)", {
  res <- count_kmers(
    rep("a", 10),
    c(1, 2),
    c("a", "b"),
    FALSE
  )
  expected <- c("a.b.a"=0, "b.b.a"=0, "b.b.b"=0, "a.a.a"=5, "b.a.b"=0, "b.a.a"=0, "a.a.b"=0, "a.b.b"=0)
  expect_mapequal(res, expected)
})

test_that("for empty alphabet empty vector should be returned", {
  expect_equal(count_kmers(c("a", "b"), c(1), c(), FALSE), c())
})

test_that("k-mers that contain letters not included in alphabet should not be generated", {
  res <- count_kmers(c("a", "b", "c", "d", "a", "b", "c"),
                     c(1,2),
                     c("a", "c", "b"),
                     TRUE)
  expect_mapequal(res, c("0_a.c.b"=1))
})
