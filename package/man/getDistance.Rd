\name{getDistance}
\alias{getDistance}

\title{Get the linear distance for each interaction}

\description{Compute the distance between interacting regions on the linear genome.}

\usage{
getDistance(data, type=c("mid", "gap", "span"))
}

\arguments{
\item{data}{a \code{DIList} object}
\item{type}{a character string specifying the type of distance to compute}
}

\value{
An integer vector is returned containing the distances between interacting regions for each pair in \code{data}.
}

\details{
The \code{getDistance} function examines each interaction in \code{data$pairs} and computes the distance between the corresponding regions on the linear genome.
Options to compute the distances are \code{mid}, for the distance between the midpoints of the regions; \code{gap}, for the distance of the interval lying between the regions (set to zero, if regions are overlapping); and \code{span}, for the total distance spanned by the interaction, including both regions. 
Interchromosomal interactions are marked with \code{NA}.
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
y <- squareCounts(fout, param, width=10, filter=-1)
getDistance(y)
getDistance(y, type="gap")
getDistance(y, type="span")

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\keyword{filtering}
