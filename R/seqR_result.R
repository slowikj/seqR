#' Bind rows of several k-mer matrices
#' 
#' @description
#' 
#' The function binds rows of several input k-mer matrices (of type \link[Matrix]{Matrix}),
#' which are results of \link[seqR]{count_kmers} and \link[seqR]{count_multimers}.
#' This implementation also handles properly k-mer matrices that do not have
#' the same columns, as opposed to the implementation of \link[base]{rbind}.
#' 
#' @param ... k-mer matrices of type \link[Matrix]{Matrix}
#' 
#' @return a k-mer matrix of type \link[Matrix]{Matrix} that is the result of the rbind operation
#' 
#' @examples
#' 
#' batch_size <- 1
#' 
#' # k-mer counting
#' resA <- count_kmers(c("AAAAA", "ASASSSSASSA"), k=5, batch_size=batch_size)
#' resB <- count_multimers(c("HWHSHS", "AASDCASD"), k_vector=c(3, 5), batch_size=batch_size)
#' 
#' # rbind
#' res <- rbind_columnwise(resA, resB)
#' 
#' @seealso Function that count k-mers of one type: count_kmers
#' @seealso Function that counts many k-mer variants in the single invocation: count_multimers
#' 
#' @md
#' @export
rbind_columnwise <- function(...) {
  input <- lapply(list(...),
                  function(x) slam::as.simple_triplet_matrix(x))
  .convert_internal_result_to_seqR_Matrix_class(.cpp_merge_kmer_results(input))
}

.convert_internal_result_to_seqR_Matrix_class <- function(seqR_list) {
  .convert_seqR_list_to_Matrix_class(seqR_list)
}

.convert_seqR_list_to_Matrix_class <- function(seqR_list) {
  if (length(seqR_list$i) == 0) {
    Matrix::Matrix(nrow=seqR_list$nrow,
                   ncol=0,
                   sparse = TRUE)
  } else {
    dimnames <- .get_dimnames(seqR_list)
    
    Matrix::sparseMatrix(
      i = seqR_list$i,
      j = seqR_list$j,
      x = seqR_list$v,
      dims = c(seqR_list$nrow, seqR_list$ncol),
      dimnames = dimnames
    )
  }
}

.get_dimnames <- function(seqR_list) {
  list(NULL, seqR_list$names)
}
