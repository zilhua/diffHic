\name{totalCounts}
\alias{totalCounts}

\title{Get the total counts}

\description{Get the total number of read pairs in a set of Hi-C libraries.}

\usage{
squareCounts(files, fragments=NULL, restrict=NULL)
}

\arguments{
	\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
	\item{fragments}{a \code{GRanges} object containing the boundaries of restriction fragments on the genome of interest}
    \item{restrict}{a character vector specifying the allowable chromosome names}
}

\value{
An integer vector is returned containing the total number of read pairs in each library.
}

\details{
As the name suggests, this function counts the total number of read pairs in each index file prepared by \code{\link{countPairs}}.
Specification of \code{fragments} is recommended to ensure that the chromosome names in each index file are consistent with those in the desired genome (e.g., from \code{\link{cutGenome}}).
Similarly, \code{restrict} can be used to restrict the chromosomes across which read pairs are counted.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, cuts, file=fout))
countPairs(fout)

# Counting totals.
totalCounts(fout)
totals(squareCounts(fout, cuts, width=10))

totalCounts(fout, restrict="chrA")
totals(squareCounts(fout, cuts, width=10, restrict="chrA"))

totalCounts(fout, restrict="chrB")
totals(squareCounts(fout, cuts, width=10, restrict="chrB"))

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{countPairs}}, 
\code{\link{cutGenome}},
}

