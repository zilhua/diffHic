\name{getDistance}
\alias{getDistance}
\alias{getArea}

\title{Get linear distance and interaction area}

\description{Compute the distance and area between interacting regions on the linear genome.}

\usage{
getDistance(data, type=c("mid", "gap", "span"))
getArea(data)
}

\arguments{
\item{data}{a named list as produced by \code{\link{squareCounts}}}
\item{type}{a character string specifying the type of distance to compute}
}

\value{
For \code{getDistance}, an integer vector is returned containing the distances between interacting regions for each pair in \code{data}.

For \code{getArea}, a numeric vector is returned containing the area in the interaction space for each pair of regions in \code{data}.
}

\details{
The \code{getDistance} function examines each interaction in \code{data$pairs} and computes the distance between the corresponding regions on the linear genome.
Options to compute the distances are \code{mid}, for the distance between the midpoints of the regions; \code{gap}, for the distance of the interval lying between the regions (set to zero, if regions are overlapping); and \code{span}, for the total distance spanned by the interaction, including both regions. 
Interchromosomal interactions are marked with \code{NA}.

The \code{getArea} function returns the area in the interaction space for each pair of regions.
This area is reported as the number of pairs of restriction fragments.
Note that any areas touching the diagonal of the interaction space are slightly overestimated.
This is because the interaction space is symmetric about the diagonal, so any information/area on the other side is redundant and should not be considered.
}

\author{Aaron Lun}

\seealso{
	\code{\link{squareCounts}}
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up the parameters
fout <- "output.h5"
preparePairs(hic.file, cuts, fout)
countPairs(fout)

# Collating to count combinations.
y<-squareCounts(fout, cuts, width=10, filter=-1)
getDistance(y)
getDistance(y, type="gap")
getDistance(y, type="span")

\dontshow{
unlink(fout, recursive=TRUE)
}
}

