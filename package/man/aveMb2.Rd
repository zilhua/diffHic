\name{aveMb2}
\alias{aveMb2}

\title{Area-adjusted average abundance}
\description{Calculate the average abundance-per-megabase squared for each bin pair.}

\usage{
aveMb2(data, areas=NULL, prior.count=2, ...) 
}

\arguments{
\item{data}{a \code{DIList} object produced by \code{\link{squareCounts}}}
\item{areas}{a numeric scalar or vector, indicating the area (in base pairs-squared) covered by each bin pair}
\item{prior.count}{an integer scalar, specifying the prior count to be used}
\item{...}{other arguments to be passed to \code{\link{aveLogCPM}}}
}

\details{
This computes the area-adjusted average abundance for each bin pair, i.e., the average abundance-per-megabase squared.
The idea is to allow comparisons between sets of bin pairs of differing widths.
For example, abundances for larger bin pairs are scaled down, while those for smaller bin pairs are scaled up.

The area for each bin pair can be specified using the \code{areas} argument.
Otherwise, it is defined using the \code{\link{getArea}} function.
Areas should be specified in terms of base pairs-squared, such that proper scaling can be performed.

Some finesse is involved in dealing with the prior count.
If the same \code{prior.count} is used for both small and large areas, downscaling of the latter will result in a smaller effective prior value being added.
This can lead to spurious differences in the final abundances.
To avoid this, the prior count is scaled up for larger areas prior to computing the average.
}

\seealso{
\code{\link{aveLogCPM}},
\code{\link{squareCounts}}
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

# Setting up the parameters (low filter is just to get more bin pairs)
fout <- "output.h5"
invisible(preparePairs(hic.file, param, fout))
y <- squareCounts(fout, param, width=50, filter=1)

aveMb2(y)
aveMb2(y, prior.count=5)
aveMb2(y, prior.count=0.1)

# Getting larger areas, to see the effect of the prior.count
areas <- getArea(y) * 1e10
aveMb2(y, areas=areas)
aveMb2(y, areas=areas, prior.count=5)
aveMb2(y, areas=areas, prior.count=0.1)
}
