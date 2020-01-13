library(testthat)

invoke_test <- function(expected_res, ...) {
  expect_equal(seqR::decode_kmer(...), expected_res)
}

test_that("decode non positional 3-mer, alfabet contains only single characters", {
  invoke_test(expected_res="a.b.c",
              encoded_sequence=c(1,2,3),
              d=c(0,0),
              begin_position=1,
              df_code2str=data.frame(code=c(1,2,3), string=c("a", "b", "c")),
              positional_kmer=FALSE)
})

test_that("decode positional 3-mer, alfabet contains only single characters", {
  invoke_test(expected_res="1_a.b.c",
              encoded_sequence=c(1,2,3),
              d=c(0,0),
              begin_position=1,
              df_code2str=data.frame(code=c(1,2,3), string=c("a", "b", "c")),
              positional_kmer=TRUE)
})

test_that("decode non positional 3-mer, alfabet contains also multiple characters", {
  invoke_test(expected_res="aa.b.cc",
              encoded_sequence=c(1,2,3),
              d=c(0,0),
              begin_position=1,
              df_code2str=data.frame(code=c(1,2,3), string=c("aa", "b", "cc")),
              positional_kmer=FALSE)
})

test_that("decode positional 3-mer, alfabet contains also multiple characters", {
  invoke_test(expected_res="1_aa.b.cc",
              encoded_sequence=c(1,2,3),
              d=c(0,0),
              begin_position=1,
              df_code2str=data.frame(code=c(1,2,3), string=c("aa", "b", "cc")),
              positional_kmer=TRUE)
})

test_that("decode positional 3-mer starting from 2, alfabet contains also multiple characters", {
  invoke_test(expected_res="2_b.cc.cc",
              encoded_sequence=c(1,2,3,3,1),
              d=c(0,0),
              begin_position=2,
              df_code2str=data.frame(code=c(1,2,3), string=c("aa", "b", "cc")),
              positional_kmer=TRUE)
})

test_that("decode positional gapped 3-mer starting from 2, alfabet contains also multiple characters", {
  invoke_test(expected_res="2_b.cc.aa",
              encoded_sequence=c(1,2,2,3,1),
              d=c(1,0),
              begin_position=2,
              df_code2str=data.frame(code=c(1,2,3), string=c("aa", "b", "cc")),
              positional_kmer=TRUE)
})
