#' Count k-mers of various types for a given collection of sequences
#' 
#' @description 
#' 
#' This is a convenient wrapper over \link[seqR]{count_kmers} function
#' in order to enable the computation of multiple types of k-mers
#' in a single invocation of the function
#' (for more information on k-mer types see \code{Supported variants of k-mers}
#' in \link[seqR]{count_kmers}).
#' 
#' A user can input multiple k-mer configurations in the following way.
#' Each parameter that is related to the configuration
#' (i.e., \code{k_vector}, \code{positional_vector}, and \code{kmer_gaps_list})
#' is represented in a sequential form (i.e., a list or a vector).
#' The i-th entry of each sequence corresponds to the i-th configuration.
#' 
#' 
#' @param sequences input sequences of one of two supported types,
#' either \code{string vector} or \code{list} of \code{string vectors}
#' 
#' @param k_vector an \code{integer vector} that represents the lengths of k-mers.
#' The i-th element corresponds to the value of \code{k} for the i-th k-mer configuration
#' 
#' @param kmer_alphabet a \code{string vector} that represents elements of a sequence
#' that should be considered
#' 
#' @param positional_vector a \code{logical vector} consisting of k-mer configurations related to the positional part.
#' The i-th element corresponds to the i-th k-mer configuration (i.e., whether the k-mer is positional or not)
#' 
#' @param kmer_gaps_list a \code{list} of \code{integer vectors} that represents the lengths of gaps of k-mers
#' for each configuration separately. The i-th element of the list corresponds to the lengths of gaps of the i-th
#' k-mer configuation
#' 
#' @param with_kmer_counts a single \code{logical} value that determines whether the result
#' should contain k-mer frequencies
#' 
#' @param with_kmer_names a single \code{logical} value that determines whether the result
#' should contain human-readable k-mer names
#' 
#' @param batch_size a single \code{integer} value representing the number of sequences
#' that are being processed in a single step
#' (for more information see section `Configurable size of batch of sequences`
#' in \link[seqR]{count_kmers})
#' 
#' @param hash_dim a single \code{integer} value representing the length of hash vector
#' that is internally used in the algorithm
#' (for more information see section `Configurable dimension of the hash value of a k-mer`
#' in \link[seqR]{count_kmers})
#' 
#' @param verbose a single \code{logical} value representing whether a user wants to get
#' extra information on the current state of computations
#' 
#' @inherit count_kmers return
#' 
#' @details
#' The comprehensive description of supported features is available
#' in the documentation of \link[seqR]{count_kmers} function.
#'  
#' @examples
#' 
#' # Counting 1-mers
#' count_multimers(c("AAAACFVV", "AAAAAA", "AAAAD"), k_vector = c(1))
#' 
#' # Counting 1-mers and 2-mers
#' count_multimers(c("AAAACFVV", "AAAAAA", "AAAAD"), k_vector = c(1, 2))
#' 
#' # Counting 1-mers, 2-mers, and gapped 2-mers with the length of the gap = 1
#' count_multimers(
#'    c("AAAACFVV", "AAAAAA", "AAAAD"),
#'    k_vector = c(1, 2)),
#'    kmer_gaps = list(NULL, NULL, c(1)))
#' 
#' # Counting 3-mers, positional 3-mers, and positional gapped 2-mers with the length of the gap = 1
#' count_multimers(
#'    c("AAAACFVV", "AAAAAA", "AAAAD"),
#'    k_vector = c(3, 3, 2),
#'    kmer_gaps_list = list(NULL, NULL, c(1)),
#'    positional_vector = c(FALSE, TRUE, TRUE))
#' 
#' @seealso \link[seqR]{count_kmers}
#' @seealso \link[seqR]{rbind_columnwise}
#' @include count_kmers.R
#' @export
count_multimers <- function(sequences,
                            k_vector,
                            kmer_alphabet = getOption("seqR_alphabet_default"),
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
      kmer_alphabet = kmer_alphabet,
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
