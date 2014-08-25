\name{boxPairs}
\alias{boxPairs}

\title{Put bin pairs into boxes}
\description{Match smaller bin pairs to the larger bin pairs in which they are nested.}

\usage{
boxPairs(reference, ..., fragments)
}

\arguments{
\item{reference}{an integer scalar specifying the reference bin size}
\item{...}{one or more named \code{DIList} objects produced by \code{\link{squareCounts}}, with smaller bin sizes than \code{reference}}
\item{fragments}{a \code{GRanges} object containing the boundaries of the restriction fragments in the genome}
}

\value{
A list object containing:
\item{indices}{a named list of integer vectors for every \code{DIList} in the ellipsis, see Details} 
\item{pairs}{a \code{DIList} describing the parent bin pairs in which the input bin pairs are nested}
}

\details{
Consider the bin size specified in \code{reference}. 
Pairs of these bins are referred to here as the parent bin pairs, and are described in the output \code{pairs} and \code{region}.
The function accepts a number of \code{DIList} objects of bin pair data in the ellipsis, referred to here as input bin pairs.
The aim is to identify the parent bin pair in which each input bin pair is nested.

One vector will be present in the output \code{indices} for each input \code{DIList} in the ellipsis.
In each vector, each entry represents an index for a single input bin pair in the corresponding \code{DIList}.
This index points to a row of \code{pairs} and indicates the parent bin pair in which that input bin pair is nested.

Note that the counts in the returned \code{pairs} are repurposed for this function.
For each parent bin pair, the count matrix in \code{pairs} specifies the number of nested bin pairs for a given input \code{DIList}.
Similarly, the \code{totals} vector contains the total number of input bin pairs in each input \code{DIList}.
All columns/entries are named according to the name of the corresponding input \code{DIList}.

All input \code{DIList} objects in the ellipsis must be constructed carefully.
In particular, the value of \code{width} in \code{\link{squareCounts}} must be such that \code{reference} is an exact multiple of each \code{width}.
This is necessary to ensure complete nesting.
Otherwise, the behaviour of the function will not be clearly defined.

Some users may wish to identify bin pairs in one \code{DIList} that are nested within bin pairs in another \code{DIList}.
This can be done by setting the larger bin size as \code{reference}, and supplying both bin pair \code{DIList} objects in the ellipsis.
It is then a trivial matter to \code{\link{match}} the output indices for the smaller bin pairs to the output indices of the larger bin pairs.
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
y2 <- squareCounts(fout, cuts, width=25, filter=1)

# Clustering.
boxPairs(reference=100, larger=y, smaller=y2, fragments=cuts)
boxPairs(reference=50, larger=y, smaller=y2, fragments=cuts)

\dontshow{
unlink(fout, recursive=TRUE)
}
}