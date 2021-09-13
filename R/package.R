#' seqR: Fast and Comprehensive K-Mer Counting Package
#' 
#' @description
#' 
#' The \code{seqR} package provides in-memory, probabilistic,
#' highly-optimized, and multi-threaded implementation of k-mer counting.
#' 
#' @author Jadwiga Słowik
#' 
#' @docType package
#' 
#' @aliases seqR
#' 
#' @name seqR-package
#' @useDynLib seqR, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @examples
#' # Load exemplary sequences
#' data(CsgA)
#' 
#' # Counting 1-mers (amino acid composition)
#' count_kmers(
#'     CsgA,
#'     k = 1,
#'     batch_size = 1)   
#' 
#' 
#' # Counting 1-mers and 2-mers
#' count_multimers(
#'     CsgA,
#'     k_vector = c(1, 2),
#'     batch_size = 1)
NULL
