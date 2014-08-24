\name{getArea}
\alias{getArea}

\title{Get interaction area}

\description{Compute area in the interaction space for each pair of regions.}

\usage{
getArea(data, fragments=NULL)
}

\arguments{
\item{data}{a \code{DIList} object}
\item{fragments}{a \code{GRanges} object containing the restriction fragment coordinates}
}

\value{
A numeric vector is returned containing the area in the interaction space for each pair of regions in \code{data}.
}

\details{
The \code{getArea} function returns the area in the interaction space for each pair of regions.
This area is reported as the number of pairs of restriction fragments.
If any partial overlaps occur between interacting regions, users should supply the restriction fragment coordinates in \code{fragments}.
This ensures that the area is computed in a manner that accounts for redundant areas due to symmetry around the diagonal.

Typically, specification of \code{fragments} is only required when \code{data} is produced by \code{\link{connectCounts}}.
This is because bins will not partially overlap in any significant manner when counts are generated with \code{\link{squareCounts}}.
}

\author{Aaron Lun}

\seealso{
	\code{\link{squareCounts}},
	\code{\link{connectCounts}}
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up the parameters
fout <- "output.h5"
preparePairs(hic.file, cuts, fout)
countPairs(fout)

# Collating to count combinations.
y<-squareCounts(fout, cuts, width=10)
getArea(y)
y<-squareCounts(fout, cuts, width=50)
getArea(y)
getArea(y, fragments=cuts)

\dontshow{
unlink(fout, recursive=TRUE)
}
}

