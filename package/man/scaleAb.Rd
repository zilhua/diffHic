\name{scaleAb}
\alias{aveMb2}

\title{Scale-adjusted average abundance}
\description{Calculate the average abundance with consideration of a specified scaling factor.}

\usage{
scaleAb(data, scale=1, prior.count=2, ...) 
}

\arguments{
\item{data}{a \code{DGEList} object, or a \code{DIList} object produced by \code{\link{squareCounts}}}
\item{scale}{a numeric scalar representing the downscaling factor of the abundances}
\item{prior.count}{an integer scalar, specifying the prior count to be used}
\item{...}{other arguments to be passed to \code{\link{aveLogCPM}}}
}

\details{
This function computes the scaled average abundance for each bin pair.
The idea is to allow comparisons between sets of bin pairs of differing widths.
For example, abundances for larger bin pairs are scaled down, while those for smaller bin pairs are scaled up.
Relative areas for the larger bin pairs should be specified as a non-unity \code{scale}.

Some finesse is involved in dealing with the prior count.
If the same \code{prior.count} is used for both small and large areas, downscaling of the latter will result in a smaller effective prior value being added.
This can lead to spurious differences in the final abundances.
To avoid this, the prior count is scaled up for larger areas prior to computing the average.
}

\value{
A numeric vector of scaled average abundances, with one entry for each bin pair.
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

small.bin <- 50
y <- squareCounts(fout, param, width=small.bin, filter=1)
large.bin <- 100
y2 <- squareCounts(fout, param, width=large.bin, filter=1)

scaleAb(y)
scaleAb(y2, scale=(large.bin/small.bin)^2)
scaleAb(y2, scale=getArea(y2)/getArea(y))

\dontshow{
unlink(fout)
}
}
