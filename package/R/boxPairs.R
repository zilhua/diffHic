boxPairs <- function(reference, ..., fragments)
# This function reports bin pairs that are nested within other bin pairs.  The
# idea is to consolidate smaller bin pairs into their larger counterparts for
# summarization of analyses involving multiple bin sizes.
#
# written by Aaron Lun
# 3 June 2014
# modified 4 August 2014
{
	if (is.null(fragments)) { stop("need restriction fragment boundaries to assemble boxes") }
	parents <- .getBinID(fragments, reference)$region
	all.hits <- list(...)
	nk <- length(all.hits)

	# Collating all results in terms of parents.
	all.a <- all.t <- all.mode <- all.idx <- list()
	num.pairs <- list()
	for (x in 1:nk) {
		current <- all.hits[[x]]
		ncur <- npairs(current)
		olap <- findOverlaps(regions(current), parents, type="within", select="first")
		if (any(is.na(olap))) { stop("smaller bins must be fully contained within larger bins") }
		
		all.a[[x]] <- olap[current@anchor.id]
		all.t[[x]] <- olap[current@target.id]
		all.mode[[x]] <- rep(x, ncur)
		all.idx[[x]] <- 1:ncur
		num.pairs[[x]] <- ncur
	}
	all.a <- unlist(all.a)
	all.t <- unlist(all.t)
	all.mode <- unlist(all.mode)
	all.idx <- unlist(all.idx)
	num.pairs <- unlist(num.pairs)

	# Ordering by anchor, target.
	o <- order(all.a, all.t)
	all.a <- all.a[o]
	all.t <- all.t[o]
	all.mode <- all.mode[o]
	all.idx <- all.idx[o]

	is.diff <- c(TRUE, diff(all.a)!=0L | diff(all.t)!=0L)
	now.index <- cumsum(is.diff)
	by.mode <- split(1:length(is.diff), all.mode)
	
	indices <- list()
	total.pairs <- sum(is.diff)
	freqs <- matrix(0L, nrow=total.pairs, ncol=nk)
	for (x in 1:nk) {
		chosen <- by.mode[[as.character(x)]]
		current.out <- integer(length(chosen))
		current.out[all.idx[chosen]] <- now.index[chosen]
		indices[[x]] <- current.out
		freqs[,x] <- tabulate(current.out, nbins=total.pairs)
	}

	# Collating all unique pairs.
	colnames(freqs) <- names(num.pairs) <- names(indices) <- names(all.hits)
	return(list(indices=indices, pairs=.DIList(counts=freqs, totals=num.pairs, 
			anchors=all.a[is.diff], targets=all.t[is.diff], regions=parents)))
}
