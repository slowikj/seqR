#' @include count_kmers.R
#' @export
count_multimers <- function(sequences,
                            k_vector,
                            alphabet = getOption("seqR_alphabet_default"),
                            positional_vector = rep(getOption("seqR_positional_default"), length(k_vector)),
                            kmer_gaps_list = rep(list(c()), length(k_vector)),
                            with_kmer_counts = getOption("seqR_with_kmer_counts_default"),
                            with_kmer_names = getOption("seqR_with_kmer_names_default"),
                            batch_size = getOption("seqR_batch_size_default"),
                            hash_dim = getOption("seqR_hash_dim_default"),
                            verbose = getOption("seqR_verbose_default")) {
  if(length(k_vector) != length(kmer_gaps_list)) {
    stop("the length of 'k_vector' must have equal length to 'kmer_gaps_list' ")
  }
  
  if(length(k_vector) != length(positional_vector)) {
    stop("the length of 'k_vector' must have equal length to 'positional vector'")
  }
  
  if(!is_bool_value(verbose)) {
    stop("verbose must be a single logical value")
  }
  
  configs_num <- length(k_vector)
  r <- do.call(cbind, lapply(1:configs_num, function(index) {
    if(verbose) {
      print(paste0("Processing sequences for ", index, " config"))
    }
    
    count_kmers(
      sequences = sequences,
      k = k_vector[index],
      alphabet = alphabet,
      positional = positional_vector[index],
      kmer_gaps = kmer_gaps_list[[index]],
      with_kmer_counts = with_kmer_counts,
      with_kmer_names = with_kmer_names,
      batch_size = batch_size,
      hash_dim = hash_dim,
      verbose = verbose)
  }))
  r
}
