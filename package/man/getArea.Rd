\name{getArea}
\alias{getArea}

\title{Get interaction area}

\description{Compute area in the interaction space for each pair of regions.}

\usage{
getArea(data, bp=TRUE)
}

\arguments{
\item{data}{a \code{DIList} object}
\item{bp}{a logical scalar indicating whether areas should be reported in base-pair terms}
}

\value{
A numeric vector is returned containing the area in the interaction space for each pair of regions in \code{data}.
}

\details{
The \code{getArea} function returns the area in the interaction space for each pair of regions.
If \code{bp=TRUE}, the area is reported in terms of squared base pairs.
This tends to be the easiest to interpret.
Otherwise, the area is reported as the number of pairs of restriction fragments.
This may be more relevant to the actual resolution of the Hi-C experiment.

Some special consideration is required for areas overlapping the diagonal.
This is because counting is only performed on one side of the diagonal, to avoid redundancy.
Base-pair areas are automatically adjusted to account for this feature, based on the presence of partial overlaps between interacting regions.

For fragment-based areas, some additional work is required to properly compute areas around the diagonal for partially overlapping regions.
This is only necessary when \code{data} is produced by \code{\link{connectCounts}}.
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
y <- squareCounts(fout, param, width=10)
getArea(y)
getArea(y, bp=FALSE)

y <- squareCounts(fout, param, width=50)
getArea(y)
getArea(y, bp=FALSE)

# Doing some extra work with arbitrary regions.
regions <- GRanges("chrA", IRanges(c(1, 100, 150), c(120, 140, 160)))
regions <- suppressWarnings(c(regions, GRanges("chrB", IRanges(50, 100))))
y2 <- connectCounts(fout, param, regions=regions, filter=1L)
getArea(y2)
getArea(y2, bp=FALSE)

\dontshow{
unlink(fout, recursive=TRUE)
}
}

