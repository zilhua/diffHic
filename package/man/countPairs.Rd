\name{countPairs}
\alias{countPairs}

\title{Count Hi-C interactions}

\description{Count the number of reads corresponding to each interaction between a pair of restriction fragments in a Hi-C library.}

\usage{
countPairs(file.in, file.out=file.in, max.length=NA, min.ingap=NA, min.outgap=NA)
}

\arguments{
	\item{file.in}{a character string specifying the path to the index file produced by \code{\link{preparePairs}}}
	\item{file.out}{a character string specifying a path to an output index file}
	\item{max.length}{an integer scalar specifying the maximum length of any sequenced DNA fragment}
	\item{min.ingap}{an integer scalar specifying the minimum distance between inward-facing reads on the same chromosome}
	\item{min.outgap}{an integer scalar specifying the minimum distance between outward-facing reads on the same chromosome}
}

\value{
An integer vector is invisibly returned, containing \code{total}, the total number of read pairs; \code{length}, the number of read pairs with fragment lengths greater than \code{max.length}; \code{ingap}, the number of inward-facing read pairs with gap distances less than \code{min.ingap}; and \code{outgap}, the number of outward-facing read pairs with gap distances less than \code{min.outgap}.

Multiple dataframe objects are also produced within the specified \code{out} file. 
Each object corresponds to a pair of chromosomes.
Each row of the dataframe corresponds to a pair of restriction fragments with at least one read pair mapped between them.
Each dataframe contains several integer fields:
\item{anchor.id}{index of the anchor restriction fragment}
\item{target.id}{index of the target restriction fragment}
\item{count}{number of read pairs spanning the interaction}
}

\details{
This function determines the number of read pairs mapped between each pair of restriction fragments in \code{file.in}.
It then stores these numbers in \code{file.out} for later reference. 
This avoids the need to recount numbers when summarization into bins is performed later, e.g., with \code{\link{squareCounts}}.

Non-\code{NA} values for \code{min.ingap} and \code{min.outgap} are designed to protect against dangling ends and self-circles, respectively. 
This is particularly true when digestion is incomplete, as said structures do not form within a single  restriction fragment and cannot be identified earlier. 
These can be removed by discarding inward- and outward-facing read pairs that are too close together.

A finite value for \code{max.length} also protects against non-specific cleavage.  
This refers to the length of the actual DNA fragment used in sequencing and is computed from the distance between each read and its nearest downstream restriction site.
Off-target cleavage will result in larger distances than expected. 

Suitable values for each parameter can be obtained with the output of \code{\link{getPairData}}. 
For example, values for \code{min.ingap} can be obtained by setting a suitable lower bound on the distribution of non-\code{NA} values for \code{gap} with \code{orientation} values of 1.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

require(rhdf5)
fout <- "temp.h5"
fout2 <- "temp2.h5"
preparePairs(hic.file, cuts, fout)
x <- countPairs(fout, fout2)
h5read(fout2, "chrA/chrA")
x <- countPairs(fout, fout2, max.length=50)
h5read(fout2, "chrA/chrA")
x <- countPairs(fout, fout2, min.ingap=50)
h5read(fout2, "chrA/chrA")
x <- countPairs(fout, fout2, min.outgap=50)
h5read(fout2, "chrA/chrA")

\dontshow{
unlink(c(fout, fout2))
}
}

\author{Aaron Lun}

\seealso{
	\code{\link{preparePairs}}, \code{\link{savePairs}}, \code{\link{getPairData}}, \code{\link{squareCounts}}
}

\references{
Jin F et al. (2013). A high-resolution map of the three-dimensional chromatin interactome in human cells. \emph{Nature} doi:10.1038/nature12644.
}
