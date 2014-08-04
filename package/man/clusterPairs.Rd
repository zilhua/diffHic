\name{clusterPairs}
\alias{clusterPairs}

\title{Cluster bin pairs}
\description{Aggregate bin pairs into local clusters for summarization.}

\usage{
clusterPairs(pairs, region, tol, upper=1e6)
}

\arguments{
\item{pairs}{a dataframe specifying the anchor/target indices of the interacting bins}
\item{region}{a \code{GRanges} object containing the region corresponding to each index}
\item{tol}{a numeric scalar specifying the maximum distance between bin pairs}
\item{upper}{a numeric scalar specifying the maximum size of each cluster}
}

\value{
An integer vector of length containing a cluster ID for each anchor/index combination in \code{pairs}.
}

\details{
Inputs for \code{pairs} and \code{region} can be taken from the output of
\code{\link{squareCounts}}. It would be prudent to filter the former based on
the average abundance of counts, in order to reduce the density of bin pairs in
the interaction space. Otherwise, clusters may not be easily interpretable.

Clustering is performed by putting a bin pair in a cluster if the smallest
Chebyshev distance to any bin pair already inside the cluster is less than
\code{tol}. This is a cross between single-linkage approaches and density-based
methods, especially after filtering removes low-density regions.

If many bin pairs are present, excessively large clusters may be formed. This
function can split the clusters such that the maximum value of any dimension of
the cluster is less than \code{upper}. This improves the spatial interpretability 
of the resulting clusters.
}

\seealso{
\code{\link{squareCounts}}
}

\author{
Aaron Lun
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up the parameters (low filter is just to get more bin pairs)
fout <- "output.h5"
preparePairs(hic.file, cuts, fout)
countPairs(fout)
y <- squareCounts(fout, cuts, width=50, filter=1)

# Clustering.
clusterPairs(y$pairs, y$region, tol=10)
clusterPairs(y$pairs, y$region, tol=50)
clusterPairs(y$pairs, y$region, tol=100)
clusterPairs(y$pairs, y$region, tol=100, upper=50)

# Repeating with a more stringent filter.
y <- squareCounts(fout, cuts, width=50, filter=2)
clusterPairs(y$pairs, y$region, tol=10)
clusterPairs(y$pairs, y$region, tol=50)
clusterPairs(y$pairs, y$region, tol=100)
clusterPairs(y$pairs, y$region, tol=100, upper=50)

unlink(fout, recursive=TRUE)
}
