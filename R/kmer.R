
#' @title Count k-mers
#' #' @export
#' 
#' @description Counts positional or non positional k-mers
#' @param seq  \code{string} matrix (each row of the matrix represents one sequence) or vector
#' @param d  \code{integer} vector representing the length of gaps between consecutive items of k-mer. See details.
#' @param alphabet  \code{string} vector representing the items which only should be used for generating k-mers
#' @param pos  \code{logical} value representing whether positional k-mers should be count
#' 
#' @return named \code{integer} vector with k-mers counts
#' 
#' @details Depending on the value of parameter \code{pos} positional or non positional k-mers can be taken
#' into account. Positional k-mers are related to its position (column index in the given \code{seq} matrix),
#' so, for example k-mer abc that starts on the 1st position is a different k-mer
#' than abc that starts on the 2nd position. 
#' As for non positional k-mers, all abc, regardless of their position, are treated the same.
#' Besides, it is important to note that during the generation of k-mers the given \code{alphabet} vector
#' is taken into account which means that we do not count k-mers that contain elements not included in the alphabet
#' and additionally, if non positional k-mers are used, all possible ones to create using the given alphabet are returned,
#' even if their count is equal to 0.
#' The format of the non positional k-mer is the following:
#' each element is separated by . char.
#' The format of the positional k-mer is similar to the above one.
#' The only thing that differs is that the prefix index_ is added to its representation of elements.
#' @examples
#' count_kmers(c("a", "b", "c"), c(0), c("a", "b"), FALSE)
count_kmers <- function(seq, d, alphabet, pos) {
  if(class(seq) != "matrix") {
    seq <- matrix(seq, nrow = 1)
  }
  alphabet <- unique(alphabet)
  
  if(length(d) == 0) {
    return (count_unigrams(seq, alphabet, pos))
  } else {
    return (count_kmers_larger_than_one(seq, d, alphabet, pos))
  } 
}

