#' @include validators.R
#' @include kmer_functions_provider.R
#' @export
count_kmers <- function(sequences,
                       k,
                       alphabet,
                       positional = FALSE,
                       kmer_gaps = c(),
                       with_kmer_counts = TRUE,
                       kmer_dictionary_name = "unordered_map",
                       batch_size = 100,
                       hash_dim = 2) {
  if (is_empty(alphabet)) {
    stop("alphabet param is empty")
  }
  
  if(is_empty(sequences)) {
    stop("sequences param is empty")
  }
  
  if(!has_integers_only(k) || k <= 0) {
    stop("k should be a positive integer")
  }
  
  if(!is.null(kmer_gaps)) {
    if(!has_integers_only(kmer_gaps)) {
      stop("gaps should be an integer vector")
    }
    if(length(kmer_gaps) >= k) {
      stop("the length of kmer_gaps vector should be at most k-1")
    }
  }
  
  if(!is_positive_integer(batch_size)) {
    stop("batch size field must be a positive integer number")
  }
  
  if(is.vector(sequences)) {
    sequences <- matrix(data=sequences, nrow=1)
  }
  
  if(!is_positive_integer(hash_dim) || hash_dim > 8) {
    stop("hash_dim is a single integer number from the range [1, 8]")
  }
  
  alphabet <- unique(alphabet)
  
  if(length(kmer_gaps) == 0) {
    invoke_contiguous_kmer_function(
      sequences=sequences,
      alphabet=alphabet,
      k=k,
      positionalKMers=positional,
      withKMerCounts=with_kmer_counts,
      kmerDictionaryName=kmer_dictionary_name,
      batchSize=batch_size,
      hashDim=hash_dim)
  } else {
    invoke_gapped_kmer_function(
      sequences=sequences,
      alphabet=alphabet,
      gaps=rep(kmer_gaps, length.out=k-1),
      positionalKMers=positional,
      withKMerCounts=with_kmer_counts,
      kmerDictionaryName=kmer_dictionary_name,
      batchSize=batch_size,
      hashDim=hash_dim)
  }
}
