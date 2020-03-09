library(testthat)

invoke_test <- function(expected_res, ...) {
  res <- seqR::count_kmers(...)
  expect_mapequal(res, expected_res)
}

test_that("count non positional 2-mers", {
  invoke_test(expected_res = c("a.a"=1, "a.b"=2, "b.a"=2),
              alphabet=c("a", "b"),
              sequence=c("a", "b", "a", "b", "a", "a"),
              k=2,
              isPositionalKMer=FALSE)
})

test_that("count positional 2-mers", {
  invoke_test(expected_res = c("1_a.b"=1, "2_b.a"=1, "3_a.b"=1, "4_b.a"=1, "5_a.a"=1),
              alphabet=c("a", "b"),
              sequence=c("a", "b", "a", "b", "a", "a"),
              k=2,
              isPositionalKMer=TRUE)
})

test_that("count non positional 1-mers", {
  invoke_test(expected_res=c("a"=3, "b"=2),
              alphabet=c("a", "b"),
              sequence=c("a", "a", "b", "a", "b"),
              k=1,
              isPositionalKMer=FALSE)
})

test_that("count positional 1-mers", {
  invoke_test(expected_res=c("1_a"=1, "2_a"=1, "3_b"=1, "4_a"=1, "5_b"=1),
              alphabet=c("a", "b"),
              sequence=c("a", "a", "b", "a", "b"),
              k=1,
              isPositionalKMer=TRUE)
})

test_that("count non positional 1-mers if some sequence items are not allowed", {
  invoke_test(expected_res=c("a"=3),
              alphabet=c("a"),
              sequence=c("a", "a", "b", "a", "b"),
              k=1,
              isPositionalKMer=FALSE)
})

test_that("count positional 1-mers if some sequence items are not allowed", {
  invoke_test(expected_res=c("1_a"=1, "2_a"=1, "4_a"=1),
              alphabet=c("a"),
              sequence=c("a", "a", "b", "a", "b"),
              k=1,
              isPositionalKMer=TRUE)
})