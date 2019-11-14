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
