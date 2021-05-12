#' @include count_kmers.R
#' @export
count_multimers <- function(sequences,
                            k_vector,
                            alphabet = "all",
                            positional_vector=rep(FALSE, length(k_vector)),
                            kmer_gaps_list=rep(list(c()), length(k_vector)),
                            with_kmer_counts = TRUE,
                            with_kmer_names = TRUE,
                            kmer_dictionary_name="unordered_map",
                            batch_size=200,
                            hash_dim=2,
                            verbose=FALSE,
                            parallel_mode=TRUE) {
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
      kmer_dictionary_name = kmer_dictionary_name,
      batch_size = batch_size,
      hash_dim = hash_dim,
      verbose = verbose,
      parallel_mode = parallel_mode)
  }))
  r
}
