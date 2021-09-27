## Test environments

* macOS-latest (release), Mac OS X 10.15.7
* windows-latest (release), Microsoft Windows Server 2019
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results

* This is a reupload.
* Fixed an error thrown by **clang-ASAN**

### NOTES

#### Possibly mis-spelled words in DESCRIPTION: Mer (3:33), mer (19:41, 19:106)

Subwords *Mer* and *mer* constitute correctly spelled parts of k-mer and k-mers words,
terms widely used in bioinformatics.

#### Installed size is 6.6Mb, sub-directories of 1Mb or more: libs 6.3Mb

The package has a larger size only if it is compiled with a g++ compiler.
In this case, the size of the compiled code is increased
due to speed performance optimizations
that take advantage of C++ templates.

#### GNU make is a SystemRequirements

The RcppParallel package, which is a dependency of the seqR package,
requires GNU make to be included in a SystemRequirements.