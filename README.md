
# seqR

<!-- badges: start -->

[![R build
status](https://github.com/slowikj/seqR/workflows/R-CMD-check/badge.svg)](https://github.com/slowikj/seqR/actions)
<!-- badges: end -->

## Overview

`seqR` is an R package for fast k-mer counting. It provides

-   highly optimized (the core algorithm is written in C++)
-   in-memory
-   probabilistic (with configurable dimensionality of a hash value used
    for storing k-mers internally),
-   multi-threaded

implementation that supports

-   various variants of k-mers (contiguous, gapped, and positional
    counterparts)
-   all biological sequences (e.g., nucleic acids and proteins)

### Exported functions

The package provides three functions that facilitate k-mer counting and
processing k-mer matrices:

-   `count_kmers` (used for counting k-mers of one type)
-   `count_multimers` (a wrapper of `count_kmers`, used for counting
    k-mers of many types)
-   `rbind_columnwise` (a helper function used for merging several k-mer
    matrices that do not have same sets of columns)

## Package installation and loading

It is possible to download the development version directly from GitHub
repository:

``` r
# install.packages("devtools")
devtools::install_github("slowikj/seqR")
```

## Supported input sequences

The package supports all types of biological sequences, in particular,
nucleic acids and proteins. The sequences can be represented as the
input in following ways:

### String vector

``` r
count_kmers(c("AADVS", "FFSFAFA"),
            k=2)
#> 2 x 9 sparse Matrix of class "dgCMatrix"
#>      A.A_0 A.D_0 D.V_0 V.S_0 S.F_0 F.S_0 F.A_0 A.F_0 F.F_0
#> [1,]     1     1     1     1     .     .     .     .     .
#> [2,]     .     .     .     .     1     1     2     1     1
```

This representation is more efficient and recommended to use. However,
it does not support multi-character alphabets, as opposed to the second
representation.

### List of string vectors

``` r
count_kmers(list(c("A", "A", "A", "AA", "V"),
                 c("VA", "A", "A", "G", "D", "D")),
            k=2)
#> 2 x 7 sparse Matrix of class "dgCMatrix"
#>      AA.V_0 A.AA_0 A.A_0 VA.A_0 G.D_0 A.G_0 D.D_0
#> [1,]      1      1     2      .     .     .     .
#> [2,]      .      .     1      1     1     1     1
```

## Supported variants of k-mers

There are four major variants of k-mers that the `seqR` package supports
(for more detail see the documentation).

### Contiguous k-mers

``` r
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            k=5)
#> 2 x 9 sparse Matrix of class "dgCMatrix"
#>      D.D.D.V.S_0.0.0.0 V.S.V.A.A_0.0.0.0 D.D.V.S.V_0.0.0.0 D.V.S.V.A_0.0.0.0
#> [1,]                 1                 1                 1                 1
#> [2,]                 .                 .                 .                 .
#>      F.S.F.S.V_0.0.0.0 F.F.S.F.S_0.0.0.0 S.F.S.V.A_0.0.0.0 S.V.A.A.A_0.0.0.0
#> [1,]                 .                 .                 .                 .
#> [2,]                 1                 1                 1                 1
#>      F.S.V.A.A_0.0.0.0
#> [1,]                 .
#> [2,]                 1
```

### Gapped k-mers

``` r
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            kmer_gaps=c(2, 0)) # gapped 3-mer template: X__XX
#> 2 x 9 sparse Matrix of class "dgCMatrix"
#>      D.V.A_2.0 D.V.S_2.0 D.S.V_2.0 V.A.A_2.0 S.A.A_2.0 F.F.S_2.0 F.A.A_2.0
#> [1,]         1         1         1         1         .         .         .
#> [2,]         .         .         .         .         1         1         1
#>      F.S.V_2.0 S.V.A_2.0
#> [1,]         .         .
#> [2,]         1         1
```

### Positional contiguous k-mers

``` r
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            k=5,
            positional=TRUE)
#> 2 x 9 sparse Matrix of class "dgCMatrix"
#>      4_V.S.V.A.A_0.0.0.0 2_D.D.V.S.V_0.0.0.0 3_D.V.S.V.A_0.0.0.0
#> [1,]                   1                   1                   1
#> [2,]                   .                   .                   .
#>      1_D.D.D.V.S_0.0.0.0 5_S.V.A.A.A_0.0.0.0 4_F.S.V.A.A_0.0.0.0
#> [1,]                   1                   .                   .
#> [2,]                   .                   1                   1
#>      2_F.S.F.S.V_0.0.0.0 1_F.F.S.F.S_0.0.0.0 3_S.F.S.V.A_0.0.0.0
#> [1,]                   .                   .                   .
#> [2,]                   1                   1                   1
```

### Positional gapped k-mers

``` r
count_kmers(c("DDDVSVAA", "FFSFSVAAA"),
            kmer_gaps=c(2, 0), # gapped 3-mer template: X__XX
            positional=TRUE)
#> 2 x 9 sparse Matrix of class "dgCMatrix"
#>      1_D.V.S_2.0 2_D.S.V_2.0 3_D.V.A_2.0 4_V.A.A_2.0 5_S.A.A_2.0 2_F.S.V_2.0
#> [1,]           1           1           1           1           .           .
#> [2,]           .           .           .           .           1           1
#>      3_S.V.A_2.0 4_F.A.A_2.0 1_F.F.S_2.0
#> [1,]           .           .           .
#> [2,]           1           1           1
```

## Result type

The result is the k-mer counting procedure is a k-mer matrix,
represented as a **sparse matrix** (see [package
Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)).

The i-th row of the result matrix corresponds to the i-th input sequence
and contains the information about particular occurrences of k-mers. On
the other hand, columns correspond to k-mers found in input sequences
and they are represented in a human-readable form (for more detail see
the documentation).

Due to the usage of a sparse matrix representation, a signification
memory optimization is achieved.

``` r
get_random_string_vector <- function(seq_num, seq_len, alphabet, seed=123) {
  set.seed(seed)
  sapply(1:seq_num,
         function(n) paste0(sample(alphabet, seq_len, replace=TRUE), collapse=""))
}

r <- count_kmers(get_random_string_vector(300, 1000, LETTERS),
                 k=4)
print(pryr::object_size(r))
#> 20.3 MB
print(pryr::object_size(as.matrix(r)))
#> 543 MB
```

Importantly, such a representation is supported by machine learning
packages such as
[ranger](https://cran.r-project.org/web/packages/ranger/index.html) and
[xgboost](https://cran.r-project.org/web/packages/xgboost/index.html).

## Extra features

In order to improve user experience and meet custom needs, several extra
features are provided.

### configurable k-mer alphabet

A k-mer alphabet specifies which elements of a sequence should be
considered during the k-mer counting procedure.

K-mer counting using the alphabet derived from input sequences:

``` r
count_kmers(c("VAAAAHDFAA", "KKKKAKAAVA"),
            k=4)
#> 2 x 14 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 14 column names 'A.A.H.D_0.0.0', 'H.D.F.A_0.0.0', 'V.A.A.A_0.0.0' ... ]]
#>                                 
#> [1,] 1 1 1 1 1 1 1 . . . . . . .
#> [2,] . . . . . . . 1 1 1 1 1 1 1
```

K-mer counting using a custom alphabet:

``` r
count_kmers(c("VAAAAHDFAA", "KKKKAKAAVA"),
            k=4,
            kmer_alphabet = c("A", "K"))
#> 2 x 6 sparse Matrix of class "dgCMatrix"
#>      A.A.A.A_0.0.0 K.K.K.K_0.0.0 K.K.A.K_0.0.0 K.A.K.A_0.0.0 K.K.K.A_0.0.0
#> [1,]             1             .             .             .             .
#> [2,]             .             1             1             1             1
#>      A.K.A.A_0.0.0
#> [1,]             .
#> [2,]             1
```

### configurable batch size

This parameter specifies how many sequences are processed in a
(multi-threaded) single step (for more detail see the documentation) and
it is available only for custom tweaking/optimization purpose. For
instance, a user can specify single-threaded mode by setting this
parameter to 1.

Below there is the comparison of speed performance between
single-threaded and multi-threaded mode.

``` r
sequences <- get_random_string_vector(seq_num = 200, seq_len = 1000, alphabet = LETTERS)
microbenchmark::microbenchmark(
  multithreaded = count_kmers(sequences, k=5, batch_size = 200),
  singlethreaded = count_kmers(sequences, k=5, batch_size = 1),
  times=11
)
#> Unit: milliseconds
#>            expr      min       lq     mean   median       uq      max neval
#>   multithreaded 162.0713 163.8967 170.9068 166.3870 167.6822 222.7021    11
#>  singlethreaded 243.9163 247.2641 248.3639 247.9051 248.2058 255.6445    11
```

### verbose mode

``` r
count_kmers(c("VAAAAHDFAA", "KKKKAKAAVA", "ASDDA"),
            k=3,
            batch_size=1,
            verbose = TRUE)
#> Start processing sequences (batch: [1-1])...
#> Start processing sequences (batch: [2-2])...
#> Start processing sequences (batch: [3-3])...
#> 3 x 17 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 17 column names 'H.D.F_0.0', 'A.H.D_0.0', 'D.F.A_0.0' ... ]]
#>                                       
#> [1,] 1 1 1 1 1 2 1 . . . . . . . . . .
#> [2,] . . . . . . . 1 1 1 1 1 2 1 . . .
#> [3,] . . . . . . . . . . . . . . 1 1 1
```

### configurable dimension of the hash value of a k-mer

As the core k-mer counting algorithm intensively takes advantage of
hashing functions, `seqR` package is a probabilistic solution. However,
due to the application of multidimensional hash values, a user can
specify the length of the hash vector, and, thus, increase the
probability of a correct result.

Nevertheless, for small values of `k` (`k <= 12`), the result is
certainly correct for the default value (i.e., `hash_size = 2`).

### compute k-mers with or without their frequencies

Compare the following two code snippets:

``` r
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_counts=FALSE)
#> 2 x 11 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 11 column names 'S.S.A_0.0', 'D.F.S_0.0', 'A.A.D_0.0' ... ]]
#>                           
#> [1,] 1 1 1 1 1 1 1 . . . .
#> [2,] . . 1 . . . . 1 1 1 1
```

``` r
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_counts=TRUE)
#> 2 x 11 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 11 column names 'S.S.A_0.0', 'D.F.S_0.0', 'A.A.D_0.0' ... ]]
#>                           
#> [1,] 1 1 1 1 1 7 1 . . . .
#> [2,] . . 1 . . . . 1 4 1 3
```

### compute a result k-mer matrix with or without human-readable k-mer (column) names

Compare the following two code snippets:

``` r
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_names=FALSE)
#> 2 x 11 sparse Matrix of class "dgCMatrix"
#>                           
#> [1,] 1 1 1 1 1 7 1 . . . .
#> [2,] . . 1 . . . . 1 4 1 3
```

``` r
count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
            k=3,
            with_kmer_names=TRUE)
#> 2 x 11 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 11 column names 'S.S.A_0.0', 'D.F.S_0.0', 'A.A.D_0.0' ... ]]
#>                           
#> [1,] 1 1 1 1 1 7 1 . . . .
#> [2,] . . 1 . . . . 1 4 1 3
```

### count k-mers of several types in a single function invocation

``` r
count_multimers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
                k_vector=c(1,2,3,4),
                verbose=TRUE)
#> [1] "Processing sequences for 1 config"
#> Start processing sequences (batch: [1-2])...
#> [1] "Processing sequences for 2 config"
#> Start processing sequences (batch: [1-2])...
#> [1] "Processing sequences for 3 config"
#> Start processing sequences (batch: [1-2])...
#> [1] "Processing sequences for 4 config"
#> Start processing sequences (batch: [1-2])...
#> 2 x 37 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 37 column names 'S', 'D', 'A' ... ]]
#>                                                                                
#> [1,] 2 1 11 1 1 1 9 1 1 1 . . 1 1 1 1 1 7 1 . . . . 1 1 1 1 1 1 5 1 . . . . . .
#> [2,] . 6  6 . . . 1 5 . . 4 1 . . 1 . . . . 1 4 1 3 . . . . . . . . 1 3 1 2 1 1
```

### merge (rbind) two k-mer matrices derived from distinct sequences

As `rbind` function does not work correctly if all input matrices have
equivalent sets of columns, `seqR` package provides a new function
(`rbind_columnwise`) to satisfy this use case:

``` r
mA <- count_kmers(c("AAAAAAADFSSAAAA", "AADADADADDAD"),
                  k=1)
mB <- count_kmers(c("VVVVVAAVA", "ADSDDD", "AAAAAV"),
                  k=1)

rbind_columnwise(mA, mB)
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>      S D  A F V
#> [1,] 2 1 11 1 .
#> [2,] . 6  6 . .
#> [3,] . .  3 . 6
#> [4,] 1 4  1 . .
#> [5,] . .  5 . 1
```

## Citation

\[TODO\]
