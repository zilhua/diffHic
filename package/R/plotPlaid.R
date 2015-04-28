plotPlaid <- function(file, param, first.region, second.region=first.region,
   	width=10000, col="red", max.count=20, xlab=NULL, ylab=NULL, 
	diag=TRUE, count=FALSE, count.args=list(), ...)
# This function takes a set of boundaries and a count directory and it generates a plaid plot of 
# the result. The plot is colour coded with heat for the desired count range (white for nothing
# and shades of 'col' according to the requested 'alpha'). Binning of restriction fragments is used
# to compute colours, using a bin width of approximately 'width'. Using overlapping rectangles
# should be avoided due to rounding imprecision when computing colours.
#
# written by Aaron Lun
# sometime in 2012.
# last modified 28 April 2015
{
	first.chr <- as.character(seqnames(first.region))
	second.chr <- as.character(seqnames(second.region))
	first.start <- start(first.region)
	first.end <- end(first.region)
	second.start <- start(second.region)
	second.end <- end(second.region)
    if (length(first.chr)!=1L) { stop("exactly one first range is required for plotting") }
    if (length(second.chr)!=1L) { stop("exactly one second range is required for plotting") }

	# Setting up the parameters
	fragments <- param$fragments
	if (!(first.chr %in% seqlevels(fragments)) || !(second.chr %in% seqlevels(fragments))) { 
		stop("anchor/target chromosome names not in cut site list") 
	}
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Setting up the boundaries.
	first.min <- max(1L, first.start)
	first.max <- min(seqlengths(fragments)[[first.chr]], first.end)
	second.min <- max(1L, second.start)
	second.max <- min(seqlengths(fragments)[[second.chr]], second.end)
	if (first.min > first.max || second.min > second.max) { stop("invalid anchor/target ranges supplied") }
	
	# Identifying the fragments in our ranges of interest (with some leeway, to ensure that edges of the plot are retained).
	keep.first <- overlapsAny(fragments, first.region, maxgap=width(first.region)/2+100)
	if (first.chr!=second.chr || first.start!=second.start || first.end!=second.end) {
		keep.second <- overlapsAny(fragments, second.region, maxgap=width(second.region)/2+100)
		keep <- keep.first | keep.second
	} else {
		keep <- keep.second <- keep.first
	}
	width <- as.integer(width) 
	new.pts <- .getBinID(fragments[keep], width)
	out.id <- integer(length(fragments))
	out.id[keep] <- new.pts$id

	# Pulling out the read pair indices from each file, and checking whether chromosome names are flipped around.
	all.dex <- .loadIndices(file, seqlevels(fragments))
	frag.by.chr <- .splitByChr(fragments)
	flipped <- FALSE

	if (!is.null(all.dex[[first.chr]][[second.chr]])) {
		current <- .baseHiCParser(TRUE, file, first.chr, second.chr, 
			chr.limits=frag.by.chr, discard=discard, cap=cap)[[1]]
	} else if (!is.null(all.dex[[second.chr]][[first.chr]])) { 
		current <- .baseHiCParser(TRUE, file, second.chr, first.chr, 
			chr.limits=frag.by.chr, discard=discard, cap=cap)[[1]]
		flipped <- TRUE
	} else { current<-data.frame(anchor.id=integer(0), target.id=integer(0)) }

	# Generating a plot.
	if (is.null(xlab)) { xlab <- first.chr }
	if (is.null(ylab)) { ylab <- second.chr }
    plot(-1, -1, xlim=c(first.min, first.max), ylim=c(second.min, second.max), xlab=xlab, ylab=ylab, type="n", ...)

	# Setting up the color function.		
	my.col<-col2rgb(col)[,1]
	colfun <- function(count) { .get.new.col(my.col, pmin(1, count/max.count)) }
	if (!nrow(current))	{ return(invisible(colfun)) }

	# Getting the counts around the area of interest, and then collating them.
   	if (flipped) {
		filter.a <- keep.second
		filter.t <- keep.first
   	} else {
 	   	filter.a <- keep.first
		filter.t <- keep.second
   	}	   
   	retain <- filter.a[current$anchor.id] & filter.t[current$target.id]
    out<-.Call(cxx_count_patch, list(current[retain,]), out.id, 1L, 
		out.id[filter.t][1L], tail(out.id[filter.t], 1L)) # First and last bin indices on the target interval.
    if (is.character(out)) { stop(out) }
	
	# Assigning coordinates to x and y-axes, while checking whether chromosome names have been flipped.
	a.ranges <- new.pts$region[out[[1]]]
	t.ranges <- new.pts$region[out[[2]]]
	if (flipped) { 
		first.ranges <- t.ranges
		second.ranges <- a.ranges
	} else {
		first.ranges <- a.ranges
		second.ranges <- t.ranges
	}

	# Summoning a function to get colours.
	all.cols <- colfun(out[[3]])
	labels <- NULL
	if (count) { labels <- out[[3]] }
	.plotDiag(first.ranges, second.ranges, all.cols, diag=diag, labels=labels, label.args=count.args)
	return(invisible(colfun))
}

###########################################################
# Helper functions.

.plotDiag <- function(xranges, yranges, colors, diag=FALSE, do.label=FALSE, labels=NULL, label.args=list()) {
	for (it in 1:2) {
		# Adding boxes (and text in those boxes, if desired).
		rect(xleft=start(xranges)-0.5, xright=end(xranges)+0.5, ybottom=start(yranges)-0.5, ytop=end(yranges)+0.5, border=NA, col=colors)
		if (!is.null(labels)) { do.call(text, c(list(x=mid(ranges(xranges)), y=mid(ranges(yranges)), labels=labels), label.args)) }

		# If we want to show elements past the diagonal for intra-chromosomal plots, we do so.
		if (!diag || as.character(seqnames(xranges[1]))!=as.character(seqnames(yranges[1]))) { break }
		offdiag <- xranges!=yranges
		if (!any(offdiag)) { break }
		temp <- xranges[offdiag]
		xranges <- yranges[offdiag]
		yranges <- temp
		colors <- colors[offdiag]
	}
	invisible(NULL)
}

.get.new.col <- function(old.colour, strength) {
	old.colour <- as.integer(old.colour)
	remnant <- 255 - old.colour
	adj <- outer(remnant, 1-strength) + old.colour
	rgb(adj[1,], adj[2,], adj[3,], maxColorValue=255)
}

###########################################################

plotDI <- function(data, fc, first.region, second.region=first.region,
	col.up="red", col.down="blue",background="grey70", zlim=NULL, 
	xlab=NULL, ylab=NULL, diag=TRUE, ...)
# This function plots differential interactions.
#
# written by Aaron Lun
# created 21 November 2014
# last modified 28 April 2015
{
	first.chr <- as.character(seqnames(first.region))
	second.chr <- as.character(seqnames(second.region))
	first.start <- start(first.region)
	first.end <- end(first.region)
	second.start <- start(second.region)
	second.end <- end(second.region)
    if (length(first.chr)!=1L) { stop("exactly one first range is required for plotting") }
    if (length(second.chr)!=1L) { stop("exactly one second range is required for plotting") }

	# Setting up the boundaries.
	first.min <- max(1L, first.start)
	first.max <- min(seqlengths(regions(data))[[first.chr]], first.end)
	second.min <- max(1L, second.start)
	second.max <- min(seqlengths(regions(data))[[second.chr]], second.end)
	if (first.min > first.max || second.min > second.max) { stop("invalid anchor/target ranges supplied") }

	# Checking that our points are consistent.
	nr <- nrow(data)
	if (nr!=length(fc)) { stop("length of fold-change vector should equal number of bin pairs") }

	# Identifying the region pairs in our ranges of interest (some stretch, to allow for partial overlaps at plot edge).
	first.keep <- overlapsAny(regions(data), first.region, maxgap=width(first.region)/2+100)
	anchor.id <- anchors(data, id=TRUE)
	target.id <- targets(data, id=TRUE)
	flipped <- FALSE

	if (first.chr!=second.chr || first.start!=second.start || first.end!=second.end) {
		second.keep <- overlapsAny(regions(data), second.region, maxgap=width(second.region)/2+100)
		keep.normal <- second.keep[target.id] & first.keep[anchor.id]
		keep.flipped <- first.keep[target.id] & second.keep[anchor.id]
		if (!any(keep.normal) && any(keep.flipped)) { flipped <- TRUE } # Checking whether chromosome names are flipped.
		keep <- keep.normal | keep.flipped # Include all bin pairs while accounting for reflection around diagonal.
	} else {
		keep <- first.keep[anchor.id] & first.keep[target.id]
	}

	if (is.null(xlab)) { xlab <- first.chr }
	if (is.null(ylab)) { ylab <- second.chr }
    plot(-1, -1, xlim=c(first.min, first.max), ylim=c(second.min, second.max), xlab=xlab, ylab=ylab, type="n", ...)
	u <- par("usr") # The coordinates of the plot area
	rect(u[1], u[3], u[2], u[4], col=background, border=NA)
	box()

	# Generating a colour-generating function.
	if (is.null(zlim)) { zlim <- max(abs(fc)) }
	up.col <- col2rgb(col.up)[,1]
	down.col <- col2rgb(col.down)[,1]
	colfun <- .forgeNewFun(zlim, up.col, down.col)
	if (!any(keep))	{ return(invisible(colfun)) }
	
	# Assigning coordinates while checking if the chromosome names are flipped.		
	current <- data[keep,]
	if (flipped) {
		first.ranges <- targets(current)
		second.ranges <- anchors(current)
	} else {
		first.ranges <- anchors(current)
		second.ranges <- targets(current)
	}	

	# Making the actual plot.
	all.cols <- colfun(fc[keep])
	.plotDiag(first.ranges, second.ranges, targets(current), all.cols, diag=diag)
	return(invisible(colfun))
}

.forgeNewFun <- function(zlim, up.col, down.col) {
	function(FC) {
		all.cols <- character(length(FC))
		pos.FC <- FC > 0				
		if (any(pos.FC)) { all.cols[pos.FC] <- .get.new.col(up.col, pmin(1, FC[pos.FC]/zlim)) }
		if (!all(pos.FC)) { all.cols[!pos.FC] <- .get.new.col(down.col, pmin(1, -FC[!pos.FC]/zlim)) }
		return(all.cols)
	}
}
