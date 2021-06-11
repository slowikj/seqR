[![R build status](https://github.com/slowikj/seqR/workflows/R-CMD-check/badge.svg)](https://github.com/slowikj/seqR/actions)

`seqR` is an R package that implements feature extraction (`k-mer counting`) of biological sequences such as proteins and nucleic acids.
It supports multiple types of sequences:
* a list of strings
* a matrix of integers
* a matrix of numbers
* a matrix of strings

Additionally, several types of k-mers are supported:
* contiguous k-mers
* gapped k-mers
* position-specific (positional) k-mers

The core algorithm is implemented in C++ and it takes advantage of parallel programming.

## Package installation and loading

```{r setup, eval=FALSE}
#install.packages('devtools') # if you do not have it installed
devtools::install_github("slowikj/seqR")
library(seqR)
```

## Usage

```{r, eval=FALSE}
count_kmers(list("aaaaacbb", "aaaaa"), 5, letters)

count_multimers(list("aaaaacb", "aaaaa"), c(3, 5), letters)
```

In order to improve memory usage, the `slam::simple_triplet_matrix` is used.
If you want to get more human-readable result, use `as.matrix` method.
