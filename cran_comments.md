## Test environments

* macOS-latest (release), Mac OS X 10.15.7
* windows-latest (release), Microsoft Windows Server 2019
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results

### NOTES

* Possibly mis-spelled words in DESCRIPTION: Mer (3:33), mer (19:41, 19:106)

Subwords *Mer* and *mer* constitute correctly spelled parts of k-mer and k-mers words,
terms widely used in bioinformatics.

* Installed size is 6.6Mb, sub-directories of 1Mb or more: libs 6.3Mb

The package has a slightly increased size only when it is compiled using g++ compiler.
The main reason is related to speed performance optimizations
that take advantage of C++ templates which increase the size of the compiled code.

* GNU make is a SystemRequirements

The RcppParallel package, which is the dependency of the seqR package,
requires GNU make to be included in a SystemRequirements.