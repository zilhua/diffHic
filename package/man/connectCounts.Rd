\name{connectCounts}
\alias{connectCounts}

\title{Count connecting read pairs}
\description{Count the number of read pairs connecting pairs of user-specified regions}

\usage{ 
connectCounts(files, fragments, regions, filter=1L, type="any")
}

\arguments{
   \item{files}{a character vector containing the paths to the count file for each library}
   \item{fragments}{a \code{GRanges} object containing the boundaries of restriction fragments on the genome of interest}
   \item{regions}{a \code{GRanges} object specifying the regions between which read pairs should be counted}
   \item{filter}{an integer scalar specifying the minimum count for each interaction}
   \item{type}{a character string specifying how restriction fragments should be assigned to regions}
}

\value{
A list is returned containing:
\item{counts}{an integer matrix containing read pair counts for each interaction between two regions}
\item{totals}{an integer vector containing the number of read pairs per library}
\item{pairs}{a dataframe of integer indices specifying the interacting regions}
\item{region}{a \code{GRanges} object of the boundaries of \code{region}, rounded to a restriction site based on \code{type}}
}

\details{
Interactions of interest are defined as those formed by pairs of elements in \code{regions}. 
The number of read pairs connecting each pair of elements can then be counted in each library. 
This can be useful for quantifying/summarizing interactions between genomic features, e.g., promoters or gene bodies.

Each row of \code{counts} contains the counts across all libraries for a particular interaction.
The interaction is defined in the corresponding row of \code{pairs}, where two indices point to the interacting regions.
Indices are denoted as anchor/target in the output \code{pairs} based on the genomic order of the corresponding \code{regions}.

The value of \code{type} feeds into \code{\link{findOverlaps}} and controls the manner in which restriction fragments are assigned to each region. 
By default, a restriction fragment is assigned to one or more regions if said fragment overlaps with any part of those regions.
This means that the boundaries of each region are expanded outwards to obtain the effective coordinates in the output \code{region}.
In contrast, setting \code{type="within"} would contract each region inwards.
}

\seealso{
\code{\link{squareCounts}},
\code{\link{findOverlaps}}
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up the parameters
fout <- "output"
preparePairs(hic.file, cuts, fout)
countPairs(fout)
regions <- GRanges("chrA", IRanges(c(1, 100, 150), c(20, 140, 160)))
regions2 <- suppressWarnings(c(regions, GRanges("chrB", IRanges(50, 100))))

# Collating to count combinations.
connectCounts(fout, cuts, regions=regions, filter=1L)
connectCounts(fout, cuts, regions=regions2, filter=1L)
connectCounts(fout, cuts, regions=regions2, filter=1L, type="within")
}

