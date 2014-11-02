countNeighbors <- function(data, flank=2) 
# This function computes the count sum for the bin pairs that are neighbouring 
# each specified bin pair in 'data'. It does so by computing the sum of all
# neighbors within 'width' on the anchor or target, and adding them. This means
# that we get a count for a 'cross' shape in the interaction space. Note that 
# for this to work, you need 'filter=1L' during count loading. It's a bit of
# a kludge but it's the simplest way to do it.
#
# written by Aaron Lun
# 23 April 2014
{
	type <- match.arg(type)
	flank <- as.integer(flank)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	rdata <- .delimitFragments(regions(data))
	last.id <- rdata$end
	first.id <- rdata$start
	names(last.id) <- names(first.id) <- rdata$chr

	np <- nrow(data)
	by.chr <- split(1:np, as.character(seqnames(anchors(data))))
	output <- matrix(0L, ncol=ncol(data), nrow=np)
   	out.n <- numeric(np)

	# Running through each pair of chromosomes.
	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, as.character(seqnames(targets(data[next.chr,]))))
		a.len <- last.id[[anchor]] - first.id[[anchor]] + 1L

		for (target in names(next.chr)) {
			current.pair <- next.chr[[target]]
			all.a <- data@anchor.id[current.pair]
			all.t <- data@target.id[current.pair]
			t.len <- last.id[[target]] - first.id[[target]] + 1L
			rel.counts <- data@counts[current.pair,,drop=FALSE]

			if (target==anchor) {
				# Doing diagonal-based background estimation.
				all.counts <- all.n <- 0L
			} else {
				# Computing neighbourhood for same anchor, different target.
				o <- order(all.a, all.t)
				collected <- .Call(cxx_collect_background, all.a[o], all.t[o], rel.counts[o,], flank, t.len)
				if (is.character(collected)) { stop(collected) }
				collected[o,] <- collected
				all.counts <- collected - rel.counts
		
				# Computing neighbourhood for same target, different anchor.
				o <- order(all.t, all.a)
				collected <- .Call(cxx_collect_background, all.t[o], all.a[o], rel.counts[o,], flank, a.len)
				if (is.character(collected)) { stop(collected) }
				collected[o,] <- collected
				all.counts <- all.counts + collected - rel.counts
				
				# Computing number of squares in the neighbourhood.
				# Code automatically adjusts on the edges to preserve the number used, so it'll
				# only change if the anchor or target chromosomes are too small.
				all.n <- min(t.len-1L, flank*2L) + min(a.len-1L, flank*2L)
			}

			output[current.pair,] <- all.counts
			out.n[current.pair] <- all.n
		}
	}

	# Returning the collected counts.
	if (any(out.n==0L)) { warning("no neighbourhoods available for some bin pairs") }
	return(list(y=DGEList(counts=output, lib.size=totals(data)), n=out.n))
}
