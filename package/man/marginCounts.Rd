\name{marginCounts}
\alias{marginCounts}

\title{Collect marginal counts for each bin}

\description{Count the number of read pairs mapped to each bin across multiple Hi-C libraries.}

\usage{
marginCounts(files, fragments, width=500000, restrict=NULL)
}

\arguments{
	\item{files}{a character vector containing paths to the index files}
	\item{fragments}{a \code{GRanges} object containing the boundaries of restriction fragments in the genome}
	\item{width}{an integer scalar specifying the width of each bin}
    \item{restrict}{a character vector specifying the allowable chromosome names}
}

\value{
A \code{DIList} object containing the marginal counts for each bin.
}

\details{
The genome is first split into non-overlapping adjacent bins of size \code{width}. 
These bins are rounded to the nearest restriction site. 
The marginal count for each bin is defined as the number of read pairs with at least one end mapped within the bin.

Even though counts do not correspond to interactions, a \code{DIList} object is still used to store the output for convenience.
Bin coordinates can be extracted as the anchor or target intervals in the output.
Each row of the output refers to a single bin in the linear genome, instead of a bin pair in the interaction space.

% Cleaner than pulling out a SummarizedExperiment from the middle of nowhere.

Larger marginal counts can be collected by increasing the \code{width} value. 
However, this comes at the cost of spatial resolution (and potentially detection power) as adjacent events in the same concatenation can no longer be distinguished. 
}

\examples{
hic.file<-system.file("exdata", "hic_sort.bam", package="diffHic")
originals<-readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up the parameters
fout <- "output.h5"
preparePairs(hic.file, originals, fout)
countPairs(fout)

# Collating to count combinations.
marginCounts(fout, originals, width=10)
marginCounts(fout, originals, width=50)
marginCounts(fout, originals, width=50, restrict="chrA")
marginCounts(fout, originals, width=100)

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\author{Aaron Lun}

\seealso{
\code{\link{squareCounts}},
\code{\link{DIList-class}}
}

