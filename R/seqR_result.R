#' @export
rbind_columnwise <- function(...) {
  input <- lapply(list(...),
                  function(x) slam::as.simple_triplet_matrix(x))
  .convert_seqR_list_to_Matrix_class(.cpp_merge_kmer_results(input))
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
