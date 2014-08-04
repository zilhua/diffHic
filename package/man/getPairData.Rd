\name{getPairData}
\alias{getPairData}

\title{Get read pair data}
\description{Extract diagnostics for each read pair from an index file}

\usage{
getPairData(file, length=TRUE, orientation=TRUE, gap=TRUE)
}

\arguments{
\item{file}{character string, specifying the path to the index file produced by \code{\link{preparePairs}}}
\item{length}{logical scalar, specifying whether the length of each sequenced DNA fragment should be extracted}
\item{orientation}{logical scalar, specifying whether the strand orientation should be extracted}
\item{gap}{logical scalar, specifying whether the insert size should be extracted}
}

\value{
A dataframe is returned containing integer fields for \code{length}, \code{orientation} and \code{gap} for each read pair 
}

\details{
This is a convenience function to extract read pair diagnostics from an index file, generated from a Hi-C library with \code{\link{preparePairs}}. 
The aim is to examine the distribution of each returned value to determine the appropriate cutoffs for \code{\link{countPairs}}. 

The \code{length} refers to the length of the DNA fragment used in sequencing.
It is computed for each read pair by adding the distance of each read to the closest restriction site in the direction of the read. 
The \code{gap} simply refers to the distance between reads on the same chromosome. 
Values for interchromosomal pairs are set to \code{NA}.

For \code{orientation}, setting 0x1 or 0x2 means that the reads designated as the ``anchor'' or ``target'' respectively are on the reverse strand.
For intrachromosomal reads, an orientation value of 1 represents inward-facing reads whereas a value of 2 represents outward-facing reads. 
}

\seealso{
\code{\link{preparePairs}},
\code{\link{countPairs}}
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

tmpf <- "gunk.h5"
preparePairs(hic.file, cuts, tmpf)
getPairData(tmpf)
getPairData(tmpf, length=FALSE)
getPairData(tmpf, gap=FALSE)
getPairData(tmpf, orientation=FALSE)

\dontshow{
unlink(tmpf)
}
}
