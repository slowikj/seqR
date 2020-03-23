library(testthat)
source("utils.R")

invoke_test <- function(expected_res, ...) {
  res <- seqR::count_gapped_kmers(...)
  expect_matrices_equal(expected_res, res)
}

test_that("test one sequence with gaps (1) not positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test(expected_res=to_matrix(c("a.a"=3, "b.b"=1, "a.b"=1)),
              alphabet=c("a", "b"),
              sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
              gaps=c(1),
              positionalKMers=FALSE)
})

test_that("test one sequence with gaps (1) positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test(expected_res=to_matrix(c("1_a.a"=1, "2_a.b"=1, "3_a.a"=1, "4_b.b"=1, "5_a.a"=1)),
              alphabet=c("a", "b"),
              sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
              gaps=c(1),
              positionalKMers=TRUE)
})

test_that("test one sequence with gapps (1,0) not positional", {
  m <- matrix(c("a", "a", "a", "b", "a", "b", "a" ),
              nrow=1)
  
  invoke_test(expected_res=to_matrix(c("b.b.a"=1, "a.a.b"=2, "a.b.a"=1)),
              alphabet=c("a", "b"),
              sequenceMatrix=to_matrix(c("a", "a", "a", "b", "a", "b", "a")),
              gaps=c(1, 0),
              positionalKMers=FALSE)
})
