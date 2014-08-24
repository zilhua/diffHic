squareCounts <- function(files, fragments, width=50000, restrict=NULL, filter=1L)
# This function collates counts across multiple experiments to get the full set of results. It takes 
# a list of lists of lists of integer matrices (i.e. a list of the outputs of convertToInteractions) and
# then compiles the counts into a list object for output. 
#
# written by Aaron Lun
{
	nlibs <- length(files)
	if (nlibs==0) { 
		stop("number of libraries must be positive")
	} else if (width < 0) { 
		stop("width must be a non-negative integer")
	} 
	if (!is.integer(width)) { width<-as.integer(width) }
	if (!is.integer(filter)) { filter<-as.integer(filter) }
	new.pts <- .getBinID(fragments, width)
	chrs <- seqlevels(fragments)
	
	# Output vectors.
	full.sizes <- integer(nlibs)
	out.counts <- list()
	out.a <- out.t <- list()
	idex <- 1L

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files)
    for (anchor in names(overall)) {
		stopifnot(anchor %in% chrs)
	    if (!is.null(restrict) && !(anchor %in% restrict)) { next }
        current <- overall[[anchor]]

		for (target in names(current)) {
			stopifnot(target %in% chrs)
	        if (!is.null(restrict) && !(target %in% restrict)) { next }

			# Extracting counts and aggregating them in C++ to obtain count combinations for each bin pair.
			pairs <- .baseHiCParser(current[[target]], files, anchor, target)
			for (lib in 1:length(pairs)) { full.sizes[lib] <- full.sizes[lib] + sum(pairs[[lib]]$count) }
            out <- .Call(cxx_count_patch, pairs, new.pts$id, filter)
			if (is.character(out)) { stop(out) }
			if (!length(out[[1]])) { next }

			# Storing counts and locations. 
			if (any(out[[1]] < out[[2]])) { stop("anchor ID should be greater than target ID") }
			out.a[[idex]] <- out[[1]]
 			out.t[[idex]] <- out[[2]]
			out.counts[[idex]] <- out[[3]]
			idex<-idex+1L
		}
	}

	# Collating all the other results.
	out.a <- unlist(out.a)
	out.t <- unlist(out.t)
	out.counts <- do.call(rbind, out.counts)
	if (!nrow(out.counts)) { 
		out.a <- out.t <- integer(0)
		out.counts <- matrix(0L, ncol=nlibs, nrow=0L)
	}

	return(.DIList(counts=out.counts, totals=full.sizes, 
		anchors=out.a, targets=out.t, regions=new.pts$region))
}

## PROOF:
# Recall the enforcement of anchor >= target. Bin pairs should technically be
# reflected around the diagonal, to ensure that all points are counted, e.g.,
# if a bin pair overlaps a point in (target, anchor) form but not in (anchor,
# target) form. However, this is not required due to the enforcement above.
# 
# Consider a point (x, y) where y > x; this implies that the target range [te,
# ts] includes 'y', and the anchor range [ae, as] includes 'x'. The anchor
# range must be above the target range (i.e., as >= ts, ae >= te). If ts <= y
# <= te and as <= x <= ae, then you can fiddle with this to obtain ts <= x <=
# te and as <= y <= ae (as x < y), i.e., the point (y, x) is also covered.

####################################################################################################

.getBinID <- function(frags, width) 
# Determines which bin each restriction fragment is in. Also records the rounded
# start and stop site for each bin. Returns a set of bin ids for each restriction
# fragment on each chromosome, as well as the coordinates of each bin.
{
	width<-as.integer(width)
	out.ids<-integer(length(frags))
	out.ranges<-list()
	last<-0L
	frag.data <- .checkFragments(frags)
	nfrags <- list() 
	
	for (x in 1:length(frag.data$chr)) {
		curindex <- frag.data$start[x]:frag.data$end[x]
		curf <- frags[curindex]
		mids <- (start(curf)+end(curf))/2
		bin.id <- as.integer((mids-0.1)/width)+1L 
		# The '-0.1' in the preceding step reduces 'mids' that are exact multiples 
		# of 'width', so each bin is from (n*width, (n+1)*width] for integer 'n'.

		stuff <- rle(bin.id)
		ns <- length(stuff$value)
		stuff$values <- 1:ns
		nfrags[[x]] <- stuff$length
		out.ids[curindex] <- inverse.rle(stuff)+last
		
		endx <- cumsum(stuff$length)
		startx <- rep(1L, ns)
		if (ns>=2L) { startx[-1] <- endx[-ns]+1L }
		out.ranges[[x]] <- GRanges(frag.data$chr[x], IRanges(start(curf[startx]), end(curf[endx])))
		last <- last+ns
	}

	# Wrapping up.
	suppressWarnings(out.ranges <- do.call(c, out.ranges))
	seqlevels(out.ranges) <- seqlevels(frags)
	seqlengths(out.ranges) <- seqlengths(frags)
	out.ranges$nfrags <- unlist(nfrags)
	return(list(id=out.ids, region=out.ranges))
}

####################################################################################################

.baseHiCParser <- function(ok, files, anchor, target)
# A courtesy function, to assist with loading counts in this function and others.
{
	overall<-list()
	for (x in 1:length(ok)) {
		if (!ok[x]) { 
			out<-data.frame(anchor.id=integer(0), target.id=integer(0), count=integer(0))
		} else {
			out <- .getPairs(files[x], anchor, target)
			check <- .Call(cxx_check_input, out$anchor.id, out$target.id, out$count)
			if (is.character(check)) { stop(check) }
		}
		overall[[x]]<-out
	}
	return(overall)
}

