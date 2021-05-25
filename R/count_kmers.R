#' @include validators.R
#' @include kmer_functions_provider.R
#' @export
count_kmers <- function(sequences,
                        k = length(kmer_gaps) + 1,
                        alphabet = getOption("seqR_alphabet_default"),
                        positional = getOption("seqR_positional_default"),
                        kmer_gaps = c(),
                        with_kmer_counts = getOption("seqR_with_kmer_counts_default"),
                        with_kmer_names = getOption("seqR_with_kmer_names_default"),
                        batch_size = getOption("seqR_batch_size_default"),
                        hash_dim = getOption("seqR_hash_dim_default"),
                        verbose = getOption("seqR_verbose_default")) {
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
  
  if(!is_positive_integer(hash_dim) || hash_dim > 8) {
    stop("hash_dim is a single integer number from the range [1, 8]")
  }
  
  if(!is_bool_value(verbose)) {
    stop("verbose must be a single logical value")
  }
  
  alphabet <- unique(alphabet)
  
  params <- rlang::env(
    alphabet=alphabet,
    positional=positional,
    with_kmer_counts=with_kmer_counts,
    with_kmer_names=with_kmer_names,
    batch_size=batch_size,
    hash_dim=hash_dim,
    verbose=verbose,
    k=k
  )
  
  if(sum(kmer_gaps) == 0) {
    .invoke_contiguous_kmer_function(
      sequences=sequences,
      params=params)
  } else {
    params$gaps = rep(kmer_gaps, length.out=k-1)
    .invoke_gapped_kmer_function(
      sequences=sequences,
      params)
  }
}
