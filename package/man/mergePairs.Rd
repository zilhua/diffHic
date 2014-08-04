\name{mergePairs}
\alias{mergePairs}

\title{Merge read pairs}

\description{Merge index files for multiple Hi-C libraries into a single output file.}

\usage{
mergePairs(files, file.out, is.count=TRUE)
}

\arguments{
	\item{files}{a character vector containing the paths to the index files to be merged}
	\item{file.out}{a character string specifying the path to the output index file}
	\item{is.count}{a logical scalar indicating whether counts should be added together}
}

\value{
A merged index file is produced at the specified location.
A \code{NULL} object is invisibly returned
}

\details{
Hi-C libraries are often split into technical replicates. 
This function facilitates the merging of said replicates into a single library for downstream processing. 
Index files listed in \code{files} should be produced by \code{\link{countPairs}} or \code{\link{preparePairs}}.
If they are generated from the latter, counts are not yet available for summation so users should set \code{is.count=FALSE}.
}

\seealso{
	\code{\link{countPairs}}, \code{\link{preparePairs}}
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <-readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

fout <- "temp.h5"
fout2 <- "temp2.h5"
fout3 <- "temp3.h5"
preparePairs(hic.file, cuts, fout)
countPairs(fout, fout2)
countPairs(fout, fout3, max.length=50)

mout <- "merged"
mergePairs(c(fout2, fout3), mout)
require(rhdf5)
h5read(fout2, "chrA/chrA")
h5read(fout3, "chrA/chrA")
h5read(mout, "chrA/chrA")

\dontshow{
unlink(c(fout, fout2, fout3), recursive=TRUE)
}
}

\author{
	Aaron Lun
}
