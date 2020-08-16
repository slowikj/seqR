[![Build Status](https://travis-ci.com/slowikj/seqR.svg?branch=master)](https://travis-ci.com/slowikj/seqR)
# seqR <img src = "man/images/logo.png" align = "right" width="120"/>

`seqR` is an R package that facilitates feature extraction (`k-mer counting`) of biological sequences such as proteins and nucleic acids.
It supports multiple types of sequences:
* `tidysq` format
* a matrix of integers
* a matrix of numbers
* a matrix of strings

Addationally, several types of k-mers are supported:
* contiguous k-mers
* gapped k-mers
* position-specific (positional) k-mers

The core algorithm is implemented in C++ and it takes advantage from parallel programming.

# Package installation and loading

```{r setup, eval=FALSE}
#install.packages('devtools') # if you do not have it installed
devtools::install_github("slowikj/seqR")
library(seqR)
```

# Usage

```{r, eval=FALSE}
count_kmers(tidysq::as.sq(c("aaaaacbb")), 5, LETTERS)

count_multimers(tidysq::as.sq(c("aaaaacb")), c(3, 5), LETTERS)
```

In order to improve memory usage, the `slam::simple_triplet_matrix` is used.
If you want to get more human-readable result, use `as.matrix` method.
