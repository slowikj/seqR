[![Build Status](https://travis-ci.org/slowikj/seqR.svg?branch=master)](https://travis-ci.org/slowikj/seqR)
[![Codecov test coverage](https://codecov.io/gh/slowikj/seqR/branch/master/graph/badge.svg)](https://codecov.io/gh/slowikj/seqR?branch=master)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/slowikj/seqR?branch=master&svg=true)](https://ci.appveyor.com/project/slowikj/seqR)
# seqR <img src = "man/images/logo.png" align = "right" width="120"/>

Fast k-mer counting is crucial in working with peptides sequences. It may take part in building features for machine learning models for peptides prediction or as a part of sequence comparison software. SeqR provides fast k-mer counting implementation in R written with Rcpp.

# Package installation and loading

Current version is available on Github, it can be build on Windows, Linux and Mac OS. For package installation from R Console the \code{devtools} package is needed. To install seqR and load it to R session the following instructions have to be executed.

```{r setup, eval=FALSE}
#install.packages('devtools') # if you do not have it installed
devtools::install_github("slowikj/seqR")
library(seqR)
```
