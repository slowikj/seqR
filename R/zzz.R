.onLoad <- function(libname, pkgname) {
  prev_options <- options()
  
  new_options <- list(
    seqR_alphabet_default = "all",
    seqR_positional_default = FALSE,
    seqR_with_kmer_counts_default = TRUE,
    seqR_with_kmer_names_default = TRUE,
    seqR_batch_size_default = 200,
    seqR_hash_dim_default = 2,
    seqR_verbose_default = FALSE
  )
  
  unset_inds <- !(names(new_options) %in% names(prev_options))
  if (any(unset_inds)) {
    options(new_options[unset_inds])
  }
  
  invisible()
}