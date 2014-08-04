\name{plotPlaid}
\alias{plotPlaid}

\title{Construct a plaid plot of interactions}

\description{Plot interactions between two sequences in a plaid format with informative colouring}

\usage{
plotPlaid(file, fragments, anchor, target=anchor, 
    width=10000, col="red", cap=20, xlab=NULL, ylab=NULL, 
	diag=TRUE, count=FALSE, count.args=list(), ...)
}

\arguments{
	\item{file}{character string specifying the path to an index file produced by \code{\link{countPairs}}}
	\item{fragments}{a \code{GRanges} object containing the boundaries of the restriction fragments on each chromosome}
	\item{anchor}{a \code{GRanges} object of length 1 specifying the anchor region}
	\item{target}{a \code{GRanges} object of length 1 specifying the target region}
	\item{width}{an integer scalar specifying the width of each bin in base pairs}
	\item{col}{any type of R color to describe the color of the plot elements}
	\item{cap}{a numeric scalar specifying the count for which the darkest color is obtained}
	\item{xlab}{character string for the x-axis label on the plot, defaults to the anchor chromosome name}
	\item{ylab}{character string for the y-axis label on the plot, defaults to the target chromosome name}
	\item{diag}{a logical scalar specifying whether boxes should be shown above the diagonal for intra-chromosomal plots}
	\item{count}{a logical scalar specifying whether the count for each bin should be plotted}
	\item{count.args}{a named list of arguments to be passed to \code{\link{text}} for plotting of bin counts}
	\item{...}{other named arguments to be passed to \code{\link{plot}}}
}	

\details{
This function constructs a plaid plot by colouring the rectangle defined by the
start/end coordinates of any pair of interacting bins within the specified
anchor/target ranges. The anchor chromosome is represented on the x-axis
whereas the target chromosome is represented on the y-axis. Base pair
coordinates are used in both cases. The square in the plotting space formed by
the pair of bins on each axis represents the interaction between those bins.

For practicality, fragments are concatenated together into a larger region of
size \code{width}. Rectangles are then defined by pairs of these regions on the
anchor/target chromosomes. The color is proportional to the number of read
pairs mapped between fragments in those regions. Larger counts can be collected
by increasing the \code{width}. This increases the range of counts (and colors)
at the cost of spatial resolution. See \code{\link{squareCounts}} for more
details. 

The resolution of colors can be controlled by varying \code{cap}. All counts
above \code{cap} will be assigned maximum opacity. Any count below \code{cap}
will be assigned a color where the transparency is proportional to the size of
the count.  Smaller values of \code{cap} will improve contrast at low counts at
the cost of contrast at higher counts. Scaling \code{cap} is recommended for
valid comparisons between libraries of different sizes (e.g., larger \code{cap}
for larger libraries).

If \code{count=TRUE}, the number of read pairs will be shown on top of each
bin. This will be slower to plot but can be useful in some cases, e.g., when
more detail is required, or when the range of colors is not sufficient to
capture the range of counts in the data. 

If \code{diag=TRUE}, boxes will also be plotted above the diagonal for
intra-chromosomal plots. This is set as the default to avoid confusion when
anchor/target ordering is misspecified, but can be turned off to reduce
redundancy in visualization around the diagonal.
}

\value{
A value of \code{NULL} is invisibly returned. A plaid plot is produced on the current graphics device. 
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
originals <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))

# Setting up parameters
fout <- "temp_saved";
preparePairs(hic.file, originals, fout)
countPairs(fout)

# Constructing intra-chromosomal plaid plots around various regions
plotPlaid(fout, originals, anchor=GRanges("chrA", IRanges(1, 100)), 
    target=GRanges("chrA", IRanges(1, 200)), width=50, diag=TRUE)
plotPlaid(fout, originals, anchor=GRanges("chrA", IRanges(1, 100)), 
    target=GRanges("chrA", IRanges(1, 200)), width=50, diag=FALSE)

# Constructing plaid plots around various regions
plotPlaid(fout, originals, anchor=GRanges("chrB", IRanges(1, 100)), 
    target=GRanges("chrA", IRanges(1, 200)), width=50)
plotPlaid(fout, originals, anchor=GRanges("chrB", IRanges(1, 100)), 
    target=GRanges("chrA", IRanges(1, 200)), width=100)

# For a hypothetical second library which is half the size of the previous one:
plotPlaid(fout, originals, anchor=GRanges("chrB", IRanges(1, 100)),
    target=GRanges("chrA", IRanges(1, 200)), width=100, cap=20, count=TRUE)
plotPlaid(fout, originals, anchor=GRanges("chrB", IRanges(1, 100)),
    target=GRanges("chrA", IRanges(1, 200)), width=100, cap=40,
    count=TRUE, count.args=list(col="blue"))

\dontshow{
unlink(fout, recursive=TRUE)
}
}

\seealso{
	\code{\link{squareCounts}}, \code{\link{countPairs}} 
}

\references{
Lieberman-Aiden E et al. (2009). Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome. \emph{Science} 326, 289-293. 
}
