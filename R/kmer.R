#' @title Amino acids
#' 
#' @return A dataframe which contains basic information about all of amino acids. 
#' Every amino acid is described by it's name, short name and letter code. 
#'
#' @export
amino_acids <- function() {
  data.frame(Name = c("Arginine", "Histidine", "Lysine",
                       "Aspartic Acid", "Glutamic Acid",
                       "Serine", "Threonine", "Aspargine",
                       "Glutamine", "Cysteine", "Selenocysteine",
                       "Glycine", "Proline", "Alanine", "Valine",
                       "Isoleucine", "Leucine", "Methionine", "Phenylalanine",
                       "Tyrosine", "Tryptophan"),
             Short = c("Arg", "His", "Lys", "Asp", "Glu", "Ser",
                       "Thr", "Asn", "Gln", "Cys", "Sec", "Gly",
                       "Pro", "Ala", "Val", "Ile", "Leu", "Met",
                       "Phe", "Tyr", "Trp"),
             Letter = c("R", "H", "K", "D", "E", "S", "T", "N",
                        "Q", "C", "U", "G", "P", "A", "V", "I",
                        "L", "M", "F", "Y", "W"))
}

#' @title Amino acids letters
#' 
#' @return A vector containing all of the valid amino acid letters
#'
#' @export
amino_letters <- function() {
  as.character(amino_acids()[, "Letter"])
}

#' @title Generate random alphabet
#'
#' @param nsize  \code{integer} number of elements in the alphabet 
#' @return Random alphabet of given size created as a subset of valid amino acid letters.
#'
#' @examples
#' generate_random_alphabet(10)
#' @export
generate_random_alphabet <- function(nsize = 10) {
  set <- amino_letters()
  set[sample(length(set), size = 10)]
}

#' @title Generate random sequence
#'
#' @param nsize  \code{integer} number of elements in the sequence
#' @return Random \code{character} object of given length created from valid amino acid letters.
#'
#' @examples
#' generate_random_seq(10)
#' @export
generate_random_seq <- function(nsize = 50) {
  set <- amino_letters()
  samp <- sample(set, nsize, replace = TRUE)
  paste(samp, collapse = '')
}

#'@title Count all mers
#'
#'@param string \code{character} with the sequence
#'@return table with number of occurrence of all mers in the sequence
#'
#'@export
#'
#'@examples
#'count_mers("AAPGAGAYY")
count_mers <- function(string) {
  l <- list()
  n <- nchar(string)
  for (j in 0:n) {
    for (i in 1:(n - j)) {
      l[j*n + i] <- substr(string, i, i+j)
    }
  }
  table(unlist(l))
}

#'@title Count kmers
#'
#'@param string \code{character} with the sequence
#'@param k \code{integer} defines length of the subsequence to be counted
#'@return table with number of occurrence of all k-mers in the sequence
#'@export
#'
#'@examples
#'count_k_mers("AAPGAGAYY", 2)
count_k_mers <- function(string, k = 2) {
  l <- list()
  n <- nchar(string)
  for (i in 1:(n - k - 1)) {
      l[i] <- substr(string, i, i+k - 1)
    }
  table(unlist(l))
}
