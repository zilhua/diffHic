\name{connectCounts}
\alias{connectCounts}

\title{Count connecting read pairs}
\description{Count the number of read pairs connecting pairs of user-specified regions}

\usage{ 
connectCounts(files, param, regions, filter=1L, type="any")
}

\arguments{
   \item{files}{a character vector containing the paths to the count file for each library}
   \item{param}{a \code{pairParam} object containing read extraction parameters}
   \item{regions}{a \code{GRanges} object specifying the regions between which read pairs should be counted}
   \item{filter}{an integer scalar specifying the minimum count for each interaction}
   \item{type}{a character string specifying how restriction fragments should be assigned to regions}
}

\value{
A \code{DIList} is returned, specifying the number of read pairs mapped between pairs of \code{regions} for each library.
}

\details{
Interactions of interest are defined as those formed by pairs of elements in \code{regions}. 
The number of read pairs connecting each pair of elements can then be counted in each library. 
This can be useful for quantifying/summarizing interactions between genomic features, e.g., promoters or gene bodies.

For a pair of intervals in \code{regions}, the interaction count is defined as the number of read pairs with one read in each interval (after rounding each interval to a fragment; see below).
To avoid reporting weak interactions, pairs can be filtered to retain only those with a count sum across all libraries above \code{filter}.
In each pair, the anchor interval is defined as that with the higher start position.
Note that the end position may not be higher, due to the possibility of nested intervals in \code{regions}.

The value of \code{type} feeds into \code{\link{findOverlaps}} and controls the manner in which restriction fragments are assigned to each region. 
By default, a restriction fragment is assigned to one or more regions if said fragment overlaps with any part of those regions.
This means that the boundaries of each region are expanded outwards to obtain the effective coordinates in the output \code{region}.
In contrast, setting \code{type="within"} would contract each region inwards.

The modified regions can be extracted from the \code{regions} slot in the output \code{DIList} object.
These will be reordered according to the new start positions.
The ordering permutation can be recovered from the \code{original} metadata field of the \code{GRanges} object.
Similarly, the number of restriction fragments assigned to each interval is stored in the \code{nfrags} metadata field.

Counting will consider the values of \code{restrict}, \code{discard} and \code{cap} in \code{param}. 
See \code{\link{pairParam}} for more details.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{findOverlaps}},
\code{\link{DIList-class}}
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

# Setting up the parameters
fout <- "output"
invisible(preparePairs(hic.file, param, fout))
regions <- c(GRanges("chrA", IRanges(c(1, 100, 150), c(20, 140, 160))), 
    GRanges("chrB", IRanges(50, 100)))

# Collating to count combinations.
con <- connectCounts(fout, param, regions=regions, filter=1L)
head(counts(con))
con <- connectCounts(fout, param, regions=regions, filter=1L, type="within")
head(counts(con))

# Still works with restriction and other parameters.
con <- connectCounts(fout, param=reform(param, restrict="chrA"), 
    regions=regions, filter=1L)
head(counts(con))
con <- connectCounts(fout, param=reform(param, discard=GRanges("chrA", IRanges(1, 50))),
    regions=regions, filter=1L)
head(counts(con))
con <- connectCounts(fout, param=reform(param, cap=1), regions=regions, filter=1L)
head(counts(con))

\dontshow{
unlink(fout)
}
}

\keyword{counting}
