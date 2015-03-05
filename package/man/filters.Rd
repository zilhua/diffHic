\name{Filtering methods}
\alias{filterDirect}
\alias{filterTrended}

\usage{
filterDirect(data, ...)
filterTrended(data, span=0.25, ...)
}

\arguments{
\item{data}{a \code{DIList} object produced by \code{\link{squareCounts}}}
\item{span}{a numeric scalar specifying the bandwidth for loess curve fitting}
\item{...}{other arguments to be passed to \code{\link{scaledAverage}}}
}

\details{
The \code{filterDirect} function implements the direct filtering strategy.
The rate of non-specific ligation is estimated as the median of average abundances from inter-chromosomal bin pairs.
This rate or some multiple thereof can be used as a minimum threshold for filtering, to keep only high-abundance bin pairs.
Some finesse is necessary when calculating the median, to consider empty parts of the interaction space that are not represented by bin pairs.

The \code{filterTrended} function implements the trended filtering strategy.
The rate of non-specific compaction is estimated by fitting a trend to the average abundances against the log-distance for all intra-chromosomal bin pairs.
This rate can then be used as a minimum threshold for filtering.
For inter-chromosomal bin pairs, the threshold is the same as that from the direct filter.

Curve fitting in \code{filterTrended} is done using \code{\link{loessFit}} with a span of \code{span}.
The bin size is also added to the distance prior to log-transformation, to avoid problems with undefined values when distances are equal to zero.

The \code{scale} argument can be passed to \code{\link{scaledAverage}}, in order to handle comparisons between bin pairs of different sizes.
Check out the user's guide for more details.
}

\value{
A list is returned containing \code{abundances}, a numeric vector with the average abundances of all bin pairs in \code{data}.
For \code{filterDirect}, the list contains a numeric scalar \code{threshold}, i.e., the non-specific ligation rate.
For \code{filterTrended}, the list contains \code{threshold}, a numeric vector containing the threshold for each bin pair; and \code{log.distances}, a numeric vector with the log-distances for each bin pair.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{scaledAverage}}
}

\author{
Aaron Lun
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

# Setting up the parameters
fout <- "output.h5"
invisible(preparePairs(hic.file, param, fout))
y <- squareCounts(fout, param, width=10, filter=-1)

# Requiring at least 3-fold change.
direct <- filterDirect(y)
keep <- direct$abundances > direct$threshold + log2(3)
y[keep,]

# Requiring to be above the threshold.
trended <- filterDirect(y)
keep <- filter$abundances > filter$threshold 
y[keep,]

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\references{
Lin, YC et al. (2012) Global changes in the nuclear positioning of genes and intra- and interdomain genomic interactions that orchestrate B cell fate. \emph{Nat. Immunol.} 13. 1196-1204
}
