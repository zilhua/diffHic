diffHic
=======

Perform a differential analysis of Hi-C data

This package provides methods for read pair processing and detection of differential interactions from Hi-C data. It is based primarily on the statistical methods in the edgeR package. Users are directed to the user's guide (found in package/inst/doc) to implement their own analyses, and to understand some of the theory behind the pipeline.

To install this package, make sure that the latest version of R is installed. Then, at the R prompt, type:
```R
source("http://bioconductor.org/biocLite.R")
useDevel()
biocLite(c('edgeR', 'Rsamtools', 'GenomicRanges', 'rhdf5'))
devtools::install_github("LTLA/diffHic/package")
```
