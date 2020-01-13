library(testthat)

invoke_test <- function(encoded_sequence, k, positional_kmer, P, M, expected_res, print_res=FALSE) {
  res <- seqR::count_kmers_hashed(encoded_sequence, k, positional_kmer, P, P^(k-1) %% M, M)
  if(print_res) {
    print(res)
  }
  expect_equal(res, expected_res)
}

test_that("test sequence with all allowed items (not positional kmers)", {
  invoke_test(encoded_sequence=c(1,2,1,2,1),
              k=2,
              positional_kmer=FALSE,
              P=101,
              M=1e9 + 9,
              expected_res = data.frame(position=c(1, 2),
                                        cnt=c(2, 2))
    )
})

test_that("test sequence with one not allowed item (not positional kmers)", {
  invoke_test(encoded_sequence = c(1,2,-1,2,1),
              k=2,
              positional_kmer=FALSE,
              P=101,
              M=1e9 + 9,
              expected_res = data.frame(position=c(1, 4),
                                        cnt=c(1, 1)))
})

test_that("test sequence with more not allowed items (not positional kmers)", {
  invoke_test(encoded_sequence = c(-1, 1, -1, 1, 2, 2, -1, -1, 1, 2, -1, 2, 1, 2, 2, -1),
              k = 3,
              positional_kmer = FALSE,
              P=101,
              M=1e9 + 9,
              expected_res = data.frame(position=c(4, 12),
                                        cnt=c(2, 1)))
})

test_that("test sequence with all allowed item (positional kmers)", {
  seq <- c(1,2,1,2,1)
  k <- 2
  last_kmer_pos <- length(seq) - k + 1
  invoke_test(encoded_sequence=seq,
              k=k,
              positional_kmer=TRUE,
              P=101,
              M=1e9 + 9,
              expected_res = data.frame(position=1:last_kmer_pos,
                                        cnt=rep(1, last_kmer_pos))
  )
})

test_that("test sequence with one not allowed item (positional kmers)", {
  invoke_test(encoded_sequence = c(1,2,-1,2,1),
              k=2,
              positional_kmer=TRUE,
              P=101,
              M=1e9 + 9,
              expected_res = data.frame(position=c(1, 4),
                                        cnt=c(1, 1)))
})

test_that("test sequence with more not allowed items (not positional kmers)", {
  invoke_test(encoded_sequence = c(-1, 1, -1, 1, 2, 2, -1, -1, 1, 2, -1, 2, 1, 2, 2, -1),
              k = 3,
              positional_kmer = TRUE,
              P=101,
              M=1e9 + 9,
              expected_res = data.frame(position=c(4, 12, 13),
                                        cnt=c(1, 1, 1)))
})
