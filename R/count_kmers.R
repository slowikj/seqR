#' Count k-mers of one, specific type for a given collection of sequences
#' 
#' @description
#' This is a in-memory, probabilistic
#' (with configurable probability of exact results,
#' for more detail see section `Configurable dimension of the hash value of a k-mer`),
#' highly-optimized, and multi-threaded implementation of the k-mer counting algorithm.
#' 
#' The function supports 
#' 1. several types of k-mers (for more information see section `Supported variants of k-mers`)
#' 2. all biological sequences (e.g., nucleic acids and proteins)
#' 3. two common in-memory representations of sequences, i.e., string vectors and list of string vectors
#' (for more information see section `Supported input sequences`)
#' 
#' Moreover, several extra features are provided
#' (for more information see corresponding `details`' subsections):
#' 1. configurable k-mer alphabet
#' (i.e., which elements of a sequence should be considered during the k-mer counting procedure)
#' 2. verbose mode
#' 3. configurable batch size (i.e., how many sequences are processed in a single step)
#' 4. configurable dimension of the hash value of a k-mer
#' 5. possibility to compute k-mers with or without their frequencies
#' 6. possibility to compute a result k-mer matrix with or without human-readable k-mer (column) names
#' 
#' @param sequences input sequences of one of two supported types,
#' either \code{string vector} or \code{list} of \code{string vectors}
#' (for more information see section `Supported input sequences`)
#' 
#' @param k an \code{integer} representing the length of a k-mer
#' 
#' @param kmer_alphabet a \code{string vector} representing the elements that should be used
#' during the construction of k-mers. By default, all elements that are present in sequences
#' are taking into account
#' 
#' @param positional a single \code{logical} value that determines whether positional k-mer
#' variant should be considered
#' (for more information on k-mer types see section `Supported variants of k-mers`)
#' 
#' @param kmer_gaps an \code{integer vector} representing the lengths of gaps between consecutive
#' k-mer elements. The length of the vector should be equal to \code{k - 1}
#' 
#' @param with_kmer_counts a single \code{logical} value that determines whether the result
#' should contain k-mer frequencies
#' 
#' @param with_kmer_names a single \code{logical} value that determines whether the result
#' should contain human-readable k-mer names
#' 
#' @param batch_size a single \code{integer} value representing the number of sequences
#' that are being processed in a single step
#' (for more information see section `Configurable size of batch of sequences`)
#' 
#' @param hash_dim a single \code{integer} value representing the length of hash vector
#' that is internally used in the algorithm
#' (for more information see section `Configurable dimension of the hash value of a k-mer`)
#' 
#' @param verbose a single \code{logical} value representing whether a user wants to get
#' extra information on the current state of computations
#' 
#' @return a \code{\link[Matrix]{Matrix}} value that represents a result k-mer space.
#' The result is a sparse matrix in order to reduce memory consumption.
#' The i-th row of the matrix represents k-mers found in the i-th input sequence.
#' Each column represents a distinct k-mer.
#' The names of columns conform to human-readable schema for k-mers,
#' if param \code{with_kmer_names = TRUE}
#' (for more information see section `Human-readable representation of k-mers`)
#' 
#' 
#' @details
#' 
#' # Supported variants of k-mers
#' 
#' A user explicitly specifies k-mer configuration using the following function params:
#' \code{k}, \code{kmer_gaps}, and \code{positional}.
#' There are four major variants of k-mers (for more information see the following subsections):
#' contiguous k-mers, gapped k-mers, positional contiguous k-mers, positional gapped k-mers.
#' 
#' ## contiguous k-mers
#' 
#' Contiguous k-mers can be defined as subwords of a fixed length.
#' 
#' A user specifies the length of the k-mer with \code{k} param.
#' The other function params should be \code{kmer_gaps = NULL} (default)
#' and \code{positional = FALSE} (default).
#' 
#' For example, 3-mers of sequence AABC are:
#' AAB and ABC.
#' 
#' ## gapped k-mers
#' 
#' Gapped k-mers can be defined as subsequences (not necessarilly contiguous) of a given sequence.
#' In particular, between two contiguous elements of the sequence there might be a gap
#' of a length specified by a user separately.
#' 
#' A user specifies the length of each gap in \code{kmer_gaps} param.
#' The other function params should be \code{k = length(kmer_gaps)} (default)
#' and \code{positional = FALSE} (default).
#' 
#' For example, gapped 3-mers with gaps' lengths 0 and 1 of sequence AABCCA are:
#' AAC, ABC, BCA.
#' 
#' ## positional contiguous k-mers
#' 
#' A positional contiguous k-mer is a subtype of a contiguous k-mer
#' with extra information about the exact (start) position of the k-mer.
#' It means that two same contiguous k-mers that start at different positions
#' of a given sequence are considered not to be equal,
#' as opposed to (non-positional) contiguous k-mers.
#' 
#' A user specifies the length of the k-mer with \code{k} param.
#' The other function params should be \code{kmer_gaps = NULL} (default)
#' and \code{positional = TRUE}.
#' 
#' For example, positional contiguous 3-mers of sequence AABCCA are:
#' 1_AAB, 2_ABC, 3_BCC, 4_CCA.
#' 
#' ## positional gapped k-mers
#' 
#' A positional gapped k-mer is a subtype of a gapped k-mer
#' with extra information about the exact (start) position of the k-mer.
#' It means that two same gapped k-mers that start at different positions
#' of a given sequence are considered not to be equal,
#' as opposed to (non-positional) gapped k-mers.
#' 
#' A user specifies the length of each gap in \code{kmer_gaps} param.
#' The other function params should be \code{k = length(kmer_gaps)} (default)
#' and \code{positional = TRUE}.
#' 
#' For example, positional gapped 3-mers with gaps' lengths 0 and 1 of sequence AABCCA are:
#' 1_AAC, 2_ABC, 3_BCA.
#' 
#' # Supported input sequences
#' 
#' The function supports all types of sequences
#' (in particular, nucleic acids and polypeptides) that are represented
#' in one of the following types:
#' 1. string vectors (e.g, \code{c("AAA", "ACDA")})
#' 2. list of string vectors (e.g., \code{list(c("A", "A", "A"), c("A", "C", "D", "A"))})
#' 
#' The first representation (\code{string vector}) is more efficient
#' and recommended as a default,
#' however, it does not support multi-character alphabets,
#' as opposed to the second representation (\code{list} of \code{string vectors}).
#' 
#' # Configurable k-mer alphabet
#' 
#' As a user might want to consider only k-mers that consist of elements
#' from a custom alphabet, they can explicitly provide the \code{vector} of elements
#' to be considered. Otherwise, all k-mers that can be derived from a given collection
#' of sequences will be computed.
#' 
#' # Verbose mode
#' 
#' Verbose mode allows a user to get extra information about the current state of computations,
#' e.g., which batch of sequences is currently processed.
#' 
#' # Configurable size of batch of sequences
#' 
#' The internal algorithm processes input sequences in consecutive batches of a given size.
#' The larger the size of a batch is, the more sequences will be used at a single step
#' in which sequences are processed concurrently.
#' Therefore, a user should be aware that in order to take full advantage of multi-threading,
#' they should set the batch size to the number that is greater or equal to the number
#' of CPU cores. However, simultaneously, this number should not be too large
#' in order to not consume too much memory, as during a single set, all sequences of a given
#' batch are encoded. Thus, they consume additional memory.
#' 
#' Generally, there are two common use cases when a user should explicitly set
#' the size of a batch size (in other cases, the default value will be appropriate):
#' 1. a user does not want to process sequences concurrently
#' (then, they should set \code{batch_size = 1})
#' 2. a user needs to tweak the memory consumption accordingly to their custom needs
#' 
#' # Configurable dimension of the hash value of a k-mer
#' 
#' As the core of the internal algorithm is based on hashing techniques,
#' the function might return approximate counts.
#' In order to significantly reduce the probability of approximate results,
#' a user can specify the dimensionality of the vector that makes up the hash value
#' of a k-mer. The longer the (hash) vector is, the more probable is that
#' it gives exact results.
#' However, a user must be aware that the more dimensions the vector has,
#' both the memory consumption and CPU time increases.
#' 
#' # Possibility to compute k-mers with or without their frequencies
#' 
#' This feature is particularly used when a user wants to get information
#' related to presence or absence of k-mers. Then, they might set this feature
#' (\code{with_kmer_counts = FALSE}) and reduce the memory consumption
#' 
#' # Possibility to compute k-mer matrix with or without human-readable names (columns)
#' 
#' The aim of this feature is to optimize memory consumption and CPU time.
#' It is particularly useful when many different k-mer models are tested
#' and there is no need to get a human-readable features.
#' 
#' # Human-readable representation of k-mers
#' 
#' Each column of the result represent a single k-mer that has the following form:
#' \deqn{[p_]s1.s2....sk_g1.g2...gk-1}
#' 
#' The \code{[p_]} value is an integer that is used only in case of positional k-mers
#' (for more information on k-mer variants see section `Supported variants of k-mers`)
#' and it indicates the exact begin position of the k-mer in a sequence.
#' The \code{s1.s2....sk} part represents consecutive elements of the k-mer.
#' Finally, the \code{g1.g2...gk-1} part represents the consecutive lengths of gaps.
#' In particular, if contiguous k-mers are considered, all elements of this part is equal to 0
#' (e.g., 0.0.0 for 4-mers). Importantly, for 1-mers, this part is not present.
#' 
#' @examples
#' 
#' # Counting 1-mers af two DNA sequences
#' count_kmers(c("ACAT", "ACC"))
#' 
#' # Counting 2-mers of two DNA sequences
#' count_kmers(c("ACAT", "ACC"), k=2)
#' 
#' # Counting positional 2-mers of two DNA sequences
#' count_kmers(c("ACAT", "ACC"), k=2, positional=TRUE)
#' 
#' # Counting positional 2-mers of two DNA sequences (second representation)
#' count_kmers(list(c("A", "C", "A", "T"), c("A", "C", "C")), k=2, positional=TRUE)
#' 
#' # Counting 3-mers of two DNA sequences, considering only A and C elements
#' count_kmers(c("ACAT", "ACC"), k=2, kmer_alphabet=c("A", "C"))
#' 
#' # Counting gapped 3-mers with lengths of gaps 1 and 2
#' count_kmers(c("ACATACTAT", "ACCCCCC"), kmer_gaps=c(1,2))
#' 
#' @seealso \link[seqR]{count_multimers}
#' @seealso \link[seqR]{rbind_columnwise}
#' @include validators.R
#' @include kmer_functions_provider.R
#' @export
#' @md
count_kmers <- function(sequences,
                        k = length(kmer_gaps) + 1,
                        kmer_alphabet = getOption("seqR_kmer_alphabet_default"),
                        positional = getOption("seqR_positional_default"),
                        kmer_gaps = c(),
                        with_kmer_counts = getOption("seqR_with_kmer_counts_default"),
                        with_kmer_names = getOption("seqR_with_kmer_names_default"),
                        batch_size = getOption("seqR_batch_size_default"),
                        hash_dim = getOption("seqR_hash_dim_default"),
                        verbose = getOption("seqR_verbose_default")) {
  if (is_empty(kmer_alphabet)) {
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
  
  kmer_alphabet <- unique(kmer_alphabet)
  
  params <- rlang::env(
    kmer_alphabet=kmer_alphabet,
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
