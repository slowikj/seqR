library(testthat)
source("utils.R")

invoke_test <- function(expectedRes, ...) {
  res <- seqR::count_kmers(kmer_dictionary_name="unordered_map",
                          batch_size=100,
                          hash_dim=2,
                          verbose=FALSE,
                          parallel_mode=TRUE,
                          ...)
  expect_matrices_equal(as.matrix(res), expectedRes)
}


nrow <- 500
ncol <- 1000000
l <- sapply(1:nrow, function(i) paste0("B", strrep("A", ncol)))

test_that("large list", {
  expectedRes <- matrix(rep(c(1, 999996), nrow),
                        nrow = nrow, byrow=TRUE)
  colnames(expectedRes) <- c( "B.A.A.A.A_0.0.0.0", "A.A.A.A.A_0.0.0.0")
  
  invoke_test(
    expectedRes = expectedRes,
    alphabet = c("A", "B"),
    sequences = l,
    k = 5,
    positional = FALSE,
    with_kmer_counts = TRUE
  )
})
