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