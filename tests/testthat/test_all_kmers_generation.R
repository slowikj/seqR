library(testthat)

invoke_test <- function(expected_res, ...) {
  res <- seqR::generate_all_kmers(...)
  expect_setequal(res, expected_res)
}

test_that("generate all 1-mers for (a,b) alphabet", {
  invoke_test(expected_res=c("a", "b"),
              k=1,
              P=101,
              M=1e9 + 9,
              decoder=c("a", "b"))
})

test_that("generate all 1-mers for (a,b,c) alphabet", {
  invoke_test(expected_res=c("a", "b", "c"),
              k=1,
              P=101,
              M=1e9 + 9,
              decoder=c("a", "b", "c"))
})

test_that("generate all 2-mers for (a,b,c) alphabet", {
  invoke_test(expected_res=c("a.a", "a.b", "a.c", "b.a", "b.b", "b.c", "c.a", "c.b", "c.c"),
              k=2,
              P=101,
              M=1e9 + 9,
              decoder=c("a", "b", "c"))
})

test_that("generate all 2-mers for (aa, bb, c) alphabet", {
  invoke_test(expected_res=c("aa.aa", "aa.bb", "aa.c", "bb.aa", "bb.bb", "bb.c", "c.aa", "c.bb", "c.c"),
              k=2,
              P=101,
              M=1e9 + 9,
              decoder=c("aa", "bb", "c"))
})

