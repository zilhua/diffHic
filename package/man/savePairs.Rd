\name{savePairs}
\alias{savePairs}

\title{Save Hi-C interactions}

\description{Save a dataframe of interactions into a directory structure for rapid chromosomal access.}

\usage{
savePairs(x, file, fragments)
}

\arguments{
	\item{x}{a sorted dataframe with integer fields \code{anchor.id} and \code{fragment.id} for each interaction}
	\item{file}{a character string specifying the path for the output index file}
    \item{fragments}{a \code{GRanges} object describing the boundaries of each restriction fragment on the reference genome}
}

\value{
Multiple Gzipped files are produced within the specified \code{dir}, as described for \code{\link{preparePairs}}.
A character string containing the path to the index file in the directory is invisibly returned.
}

\details{
This function facilitates the input of processed Hi-C data from other sources into the current pipeline. 
Each entry in \code{x} refers to the index of \code{fragments} to denote the interacting regions. 
Any additional fields in \code{x} will also be saved to file, e.g., \code{gap}, \code{length}, \code{orientation}.

The \code{x} object will be resorted by \code{anchor.id}, then \code{target.id}.
If necessary, anchor and target IDs will be switched such that the former is never less than the latter. 
Also, the coordinates of the restriction fragment boundaries in \code{fragments} must correspond to the reference genome being used. 
In most cases, these can be generated using the \code{\link{cutGenome}} function from any given \code{BSgenome} object.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <-readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

n <- 1000
all.a <- as.integer(runif(n, 1L, length(cuts)))
all.t <- as.integer(runif(n, 1L, length(cuts)))
x <- data.frame(anchor.id=pmax(all.a, all.t), target.id=pmin(all.a, all.t),
    orientation=round(runif(n, 0, 4)), length=as.integer(runif(n, 100, 10000)))
x <- x[order(x$anchor.id, x$target.id),]

fout <- "temp2.h5"
savePairs(x, fout, cuts)
require(rhdf5)
head(h5read(fout, "chrA/chrA"))

\dontshow{
unlink(fout)
}
}

\author{Aaron Lun}

\seealso{
	\code{\link{preparePairs}}, \code{\link{cutGenome}}
}

