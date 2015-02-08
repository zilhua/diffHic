\name{totalCounts}
\alias{totalCounts}

\title{Get the total counts}

\description{Get the total number of read pairs in a set of Hi-C libraries.}

\usage{
totalCounts(files, param, cap=NA)
}

\arguments{
	\item{files}{a character vector containing paths to the index files generated from each Hi-C library}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
    \item{cap}{an integer scalar or vector, indicating the maximum number of read pairs to consider for each restriction fragment pair}
}

\value{
An integer vector is returned containing the total number of read pairs in each library.
}

\details{
As the name suggests, this function counts the total number of read pairs in each index file prepared by \code{\link{preparePairs}}.
Use of \code{param$fragments} ensures that the chromosome names in each index file are consistent with those in the desired genome (e.g., from \code{\link{cutGenome}}).
Similarly, \code{param$restrict}, \code{param$discard} and \code{cap} will determine which and how many read pairs are used.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, file=fout))

# Counting totals, and comparing them.
totalCounts(fout, param)
squareCounts(fout, param, width=10)$totals

new.param <- reform(param, restrict="chrA")
totalCounts(fout, new.param)
squareCounts(fout, new.param, width=10)$totals

new.param <- reform(param, discard=GRanges("chrA", IRanges(1, 50)))
totalCounts(fout, new.param)
squareCounts(fout, new.param, width=10)$totals

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{preparePairs}}, 
\code{\link{cutGenome}},
\code{\link{pairParam}},
\code{\link{squareCounts}}
}

