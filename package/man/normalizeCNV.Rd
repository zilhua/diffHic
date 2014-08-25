\name{normalizeCNV}
\alias{normalizeCNV}

\title{Normalize CNV biases}
\description{Compute normalization offsets to remove CNV-driven and abundance-dependent biases}

\usage{
normalizeCNV(data, margins, ref.col=1, prior.count=3, split=TRUE, abundance=TRUE, ...)
}

\arguments{
\item{data}{a \code{DIList} object produced by \code{\link{squareCounts}}}
\item{margins}{a \code{DIList} object produced by \code{\link{marginCounts}}}
\item{ref.col}{an integer scalar, specifying the library to use as the reference}
\item{prior.count}{a numeric scalar specifying the prior count to use in computing marginal log-ratios}
\item{split}{a logical scalar indicating whether the marginal log-ratios should be kept as separate covariates}
\item{abundance}{a logical scalar indicating whether the average abundance should be used as a covariate}
\item{...}{other arguments to pass to \code{\link{locfit}}}
}

\details{
Each bin pair in \code{data} is associated with three covariates.
The first two are the marginal log-ratios of the corresponding bins, i.e., the log-ratio of the marginal counts between two libraries.
These represent the relative CNVs in the interacting regions between libraries.
The first covariate is the larger marginal log-ratio, whereas the second is the smaller.
The third covariate is the average abundance across all libraries.

Each bin pair is also associated with a response, i.e., the log-ratio of the interaction counts between two libraries.
A loess-like surface is fitted to the response against the three covariates, using the \code{\link{locfit}} function.
The aim is to eliminate systematic differences between libraries at any combination of covariate values.
In addition, the use of the average abundance simultaneously removes trended biases.
The fitted value can then be used as a GLM offset for each bin pair.

The \code{DIList} objects in \code{data} and \code{margins} should be constructed with the same parameters in their respective functions.
This ensures that the regions are the same, so that the marginal counts can be directly used.
Note that the marginal counts are not directly computed from \code{data} as filtering of bin pairs may be performed beforehand.

Normalization offsets are computed between libraries in a pairwise manner.
For datasets with more than two libraries, one library is arbitrarily chosen as the reference library according to \code{ref.col}.
All other libraries are then normalized to this reference library.

If \code{split=FALSE}, the two marginal log-ratios for each bin pair will be added into a single covariate.
This improves stability by reducing the covariate space, though it assumes that such addition is sensible, i.e., different combinations of log-ratios with the same sum have the same effect.
If \code{abundance=FALSE}, the average abundance will not be used as a covariate.
This also reduces the covariate space and can be used when trended biases are not present.
}

\value{
A numeric matrix of the same dimensions as \code{counts(data)}, containing log-based GLM offsets. 
}

\seealso{
\code{\link{locfit}},
\code{\link{squareCounts}},
\code{\link{marginCounts}}
}

\author{
Aaron Lun
}

\examples{
# Dummying up some data.
set.seed(3423746)
npts <- 100
npairs <- 5000
nlibs <- 4
anchors <- sample(npts, npairs, replace=TRUE)
targets <- sample(npts, npairs, replace=TRUE)
data <- diffHic:::.DIList(counts=matrix(rpois(npairs*nlibs, runif(npairs, 10, 100)), nrow=npairs),
	totals=runif(nlibs, 1e6, 2e6), anchors=pmax(anchors, targets), targets=pmin(anchors, targets),
	regions=GRanges("chrA", IRanges(1:npts, 1:npts)))
margins <- diffHic:::.DIList(counts=matrix(rpois(npts*nlibs, 100), nrow=npts),
	totals=totals(data), anchors=1:npts, targets=1:npts, regions=regions(data))

# Running normalizeCNV.
head(normalizeCNV(data, margins))
head(normalizeCNV(data, margins, prior.count=1))
head(normalizeCNV(data, margins, ref.col=3))
head(normalizeCNV(data, margins, split=FALSE))
head(normalizeCNV(data, margins, abundance=FALSE))

# Occasionally locfit will complain; increase maxk to compensate.
data@counts <- matrix(as.integer(rpois(npairs*nlibs, 20)), nrow=npairs)
tryCatch(head(normalizeCNV(data, margins)), error=function(e) e)
head(normalizeCNV(data, margins, maxk=1000))
}
