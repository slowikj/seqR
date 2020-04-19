#' @include validators.R
#' @include kmer_functions_provider.R
#' @export
count_kmers <- function(sequences,
                        k,
                        alphabet,
                        positional = FALSE,
                        kmer_gaps = c()) {
  if (is_empty(alphabet)) {
    stop("alphabet param is empty")
  }
  
  if(is_empty(sequences)) {
    stop("sequences param is empty")
  }
  
  if(is.vector(sequences)) {
    sequences <- matrix(data=sequences, nrow=1)
  }
  
  alphabet <- unique(alphabet)
  
  if(length(kmer_gaps) == 0) {
    invoke_contiguous_kmer_function(
      sequences=sequences, alphabet=alphabet, k=k, positionalKMers=positional)
  } else {
    invoke_gapped_kmer_function(
      sequences=sequences, alphabet=alphabet, gaps=rep(kmer_gaps, length.out=k-1), positionalKMers=positional)
  }
}
