
# seqR - fast and comprehensive k-mer counting package

<!-- badges: start -->

[![R build
status](https://github.com/slowikj/seqR/workflows/R-CMD-check/badge.svg)](https://github.com/slowikj/seqR/actions)
<!-- badges: end -->

## About

`seqR` is an R package for fast k-mer counting. It provides

-   **highly optimized** (the core algorithm is written in C++)
-   **in-memory**
-   **probabilistic** (with configurable dimensionality of a hash value
    used for storing k-mers internally),
-   **multi-threaded**

implementation that supports

-   **various variants of k-mers** (contiguous, gapped, and positional
    counterparts)
-   **all biological sequences** (e.g., nucleic acids and proteins)

Moreover, the result optimizes memory consumption by the application of
**sparse matrices**, compatible with machine learning packages such as
[ranger](https://cran.r-project.org/web/packages/ranger/index.html) and
[xgboost](https://cran.r-project.org/web/packages/xgboost/index.html).

## How to…

### How to install

It is possible to download the development version directly from GitHub
repository:

``` r
# install.packages("devtools")
devtools::install_github("slowikj/seqR")
```

### How to use

The package provides two functions that facilitate k-mer counting

-   `count_kmers` (used for counting k-mers of one type)
-   `count_multimers` (a wrapper of `count_kmers`, used for counting
    k-mers of many types in a single invocation of the function)

and one function used for custom processing of k-mer matrices:

-   `rbind_columnwise` (a helper function used for merging several k-mer
    matrices that do not have same sets of columns)

To learn more, see features overview vignette and documentation.

#### Examples

##### counting 5-mers

``` r
count_kmers(sequences=c("AAAAAVVAVFF", "DFGSADFGSA"),
            k=5)
#> 2 x 12 sparse Matrix of class "dgCMatrix"
#>    [[ suppressing 12 column names 'A.A.V.V.A_0.0.0.0', 'A.A.A.A.A_0.0.0.0', 'V.A.V.F.F_0.0.0.0' ... ]]
#>                             
#> [1,] 1 1 1 1 1 1 1 . . . . .
#> [2,] . . . . . . . 1 1 1 2 1
```

##### counting gapped 5-mers with gaps (0, 1, 0, 2) (XX\_XX\_\_X)

``` r
count_kmers(sequences=c("AAAAAVVAVFF", "DFGSADFGSA"),
            kmer_gaps=c(0, 1, 0, 2))
#> 2 x 7 sparse Matrix of class "dgCMatrix"
#>      A.A.V.A.F_0.1.0.2 A.A.A.A.A_0.1.0.2 A.A.V.V.F_0.1.0.2 A.A.A.V.V_0.1.0.2
#> [1,]                 1                 1                 1                 1
#> [2,]                 .                 .                 .                 .
#>      F.G.A.D.S_0.1.0.2 G.S.D.F.A_0.1.0.2 D.F.S.A.G_0.1.0.2
#> [1,]                 .                 .                 .
#> [2,]                 1                 1                 1
```

### How to cite

\[TODO\]

## Benchmarks

\[TODO\]
