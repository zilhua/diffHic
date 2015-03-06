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
# Setting up the object.
a <- 10
b <- 20
regions <- GRanges(rep(c("chrA", "chrB"), c(a, b)), IRanges(c(1:a, 1:b), c(1:a, 1:b)))

set.seed(3423)
all.anchors <- sample(length(regions), 50, replace=TRUE)
all.targets <- as.integer(runif(50, 1, all.anchors+1))
y <- DIList(matrix(0, 50, 1), anchors=all.anchors, targets=all.targets, 
    regions=regions, exptData=List(width=1))

# Clustering; note, small tolerances are used in this toy example.
clusterPairs(y, tol=1)
clusterPairs(y, tol=3)
clusterPairs(y, tol=5)
clusterPairs(y, tol=5, upper=5)
}

\keyword{clustering}
