\name{getArea}
\alias{getArea}

\title{Get interaction area}

\description{Compute area in the interaction space for each pair of regions.}

\usage{
getArea(data, param=NULL)
}

\arguments{
\item{data}{a \code{DIList} object}
\item{param}{a \code{pairParam} object containing read extraction parameters}
}

\value{
A numeric vector is returned containing the area in the interaction space for each pair of regions in \code{data}.
}

\details{
The \code{getArea} function returns the area in the interaction space for each pair of regions.
This area is reported as the number of pairs of restriction fragments.
If any partial overlaps occur between interacting regions, users should supply the restriction fragment coordinates in \code{param$fragments}.
This ensures that the area is computed in a manner that accounts for redundant areas due to symmetry around the diagonal.

Typically, specification of \code{param} is only required when \code{data} is produced by \code{\link{connectCounts}}.
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
param <- pairParam(cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, fout))

# Collating to count combinations.
y<-squareCounts(fout, param, width=10)
getArea(y, param)
y<-squareCounts(fout, param, width=50)
getArea(y)
getArea(y, param)

\dontshow{
unlink(fout, recursive=TRUE)
}
}

