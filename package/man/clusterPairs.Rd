\name{clusterPairs}
\alias{clusterPairs}

\title{Cluster bin pairs}
\description{Aggregate bin pairs into local clusters for summarization.}

\usage{
clusterPairs(data, tol, upper=1e6)
}

\arguments{
\item{data}{a \code{DIList} object}
\item{tol}{a numeric scalar specifying the maximum distance between bin pairs}
\item{upper}{a numeric scalar specifying the maximum size of each cluster}
}

\value{
An integer vector of length containing a cluster ID for each anchor/index combination in \code{pairs}.
}

\details{
Clustering is performed by putting a bin pair in a cluster if the smallest Chebyshev distance to any bin pair already inside the cluster is less than \code{tol}. 
This is a cross between single-linkage approaches and density-based methods, especially after filtering removes low-density regions.
In this manner, adjacent bin pairs in the interaction space can be clustered together.

The input \code{data} object can be taken from the output of \code{\link{squareCounts}}. 
It would be prudent to filter the former based on the average abundances, to reduce the density of bin pairs in the interaction space. 
Otherwise, clusters may be too large to be easily interpreted.

Alternatively, to avoid excessively large clusters, this function can also split each cluster into roughly-equally sized subclusters.
The maximum value of any dimension of the subclusters is approxiamtely equal to \code{upper}. 
This aims to improve the spatial interpretability of the clustering result.

% Having evenly sized clusters provides best spatial resolution and interpretability, e.g., much like binary searching.
% However, it's not clear how you can provide an even size while enforcing a strict upper bound of `upper`.
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
param <- pairParam(cuts)

# Setting up the parameters (low filter is just to get more bin pairs)
fout <- "output.h5"
invisible(preparePairs(hic.file, param, fout))
y <- squareCounts(fout, param, width=50, filter=1)

# Clustering.
clusterPairs(y, tol=10)
clusterPairs(y, tol=50)
clusterPairs(y, tol=100)
clusterPairs(y, tol=100, upper=50)

# Repeating with a more stringent filter.
y <- squareCounts(fout, param, width=50)
clusterPairs(y, tol=10)
clusterPairs(y, tol=50)
clusterPairs(y, tol=100)
clusterPairs(y, tol=100, upper=50)

\dontshow{
unlink(fout, recursive=TRUE)
}
}
