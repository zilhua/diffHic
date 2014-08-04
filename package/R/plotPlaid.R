plotPlaid <- function(file, fragments, anchor, target=anchor, 
   	width=10000, col="red", cap=20, xlab=NULL, ylab=NULL, 
	diag=TRUE, count=FALSE, count.args=list(), ...)
# This function takes a set of boundaries and a count directory and it generates a plaid plot of 
# the result. The plot is colour coded with heat for the desired count range (white for nothing
# and shades of 'col' according to the requested 'alpha'). Binning of restriction fragments is used
# to compute colours, using a bin width of approximately 'width'. Using overlapping rectangles
# should be avoided due to rounding imprecision when computing colours.
#
# written by Aaron Lun
# sometime in 2012.
{
	if (!is.integer(width)) { width<-as.integer(width) }
	achr <- as.character(seqnames(anchor))
	tchr <- as.character(seqnames(target))
	astart <- start(anchor)
	aend <- end(anchor)
	tstart <- start(target)
	tend <- end(target)
    if (length(achr)!=1L) { stop("exactly one anchor range is required for plotting") }
    if (length(tchr)!=1L) { stop("exactly one target range is required for plotting") }
	if (!(achr %in% seqlevels(fragments)) || !(tchr %in% seqlevels(fragments))) { stop("anchor/target chromosome names not in cut site list") }

	# Setting up the boundaries.
	a.min <- max(1L, astart)
	a.max <- min(seqlengths(fragments)[[achr]], aend)
	t.min <- max(1L, tstart)
	t.max <- min(seqlengths(fragments)[[tchr]], tend)
	if (a.min >= a.max || t.min >= t.max) { stop("invalid anchor/target ranges supplied") }
	
	# Identifying the fragments in our ranges of interest (with some leeway, to ensure that edges of the plot are retained).
	keep.a <- keep.t <- overlapsAny(fragments, anchor, maxgap=width(anchor)/2)
	if (anchor!=target) { 
		keep.t <- overlapsAny(fragments, target, maxgap=width(target)/2)
	}
	keep <- keep.a | keep.t
	new.pts <- .getBinID(fragments[keep], width)
	out.id <- integer(length(fragments))
	out.id[keep] <- new.pts$id

	# Pulling out the read pair indices from each file.
	all.dex <- .loadIndices(file)
	flipped <- FALSE
	if (!is.null(all.dex[[achr]][[tchr]])) {
		current <- .getPairs(file, achr, tchr)
	} else if (!is.null(all.dex[[tchr]][[achr]])) { 
		current <- .getPairs(file, tchr, achr)
		flipped <- TRUE
	} else { current<-data.frame(anchor.id=integer(0), target.id=integer(0), count=integer(0)) }

	# Generating a plot.
	if (is.null(xlab)) { xlab <- achr }
	if (is.null(ylab)) { ylab <- tchr }
	plot(-1, -1, xlim=c(a.min, a.max), ylim=c(t.min, t.max), xlab=xlab, ylab=ylab, type="n", bg="transparent", ...)
	if (!nrow(current))	{ next }

	# Getting the counts around the area of interest, and then collating them.
   	if (flipped) {
		filter.a <- keep.t
		filter.t <- keep.a
   	} else {
 	   	filter.a <- keep.a
		filter.t <- keep.t	
   	}	   
   	retain <- filter.a[current$anchor.id] & filter.t[current$target.id]
    out<-.Call(cxx_count_patch, list(current[retain,]), out.id, 1L)
    if (is.character(out)) { stop(out) }
	
	# Plotting each bin (with or without the count text).
	if (flipped) { 
		targets <- new.pts$region[out[[1]]]
		anchors <- new.pts$region[out[[2]]]
	} else {
		anchors <- new.pts$region[out[[1]]]
		targets <- new.pts$region[out[[2]]]
	}
	cur.counts<-out[[3]]
	my.col<-col2rgb(col)[,1]
	for (it in 1:2) {
		rect(xleft=start(anchors), xright=end(anchors), ybottom=start(targets), ytop=end(targets), border=NA, 
			col=rgb(my.col[1], my.col[2], my.col[3], alpha=255*pmin(1, cur.counts/cap), maxColorValue=255))
		if (count) { do.call(text, c(list(x=mid(ranges(anchors)), y=mid(ranges(targets)), labels=cur.counts), count.args)) }

		# If we want to show elements past the diagonal for intra-chromosomal plots, we do so.
		if (!diag || achr!=tchr) { break }
		offdiag <- anchors!=targets
		if (!any(offdiag)) { break }
		temp <- anchors[offdiag]
		anchors <- targets[offdiag]
		targets <- temp
		cur.counts <- cur.counts[offdiag]
	}
	return(invisible(NULL));
}

