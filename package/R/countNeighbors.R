countNeighbors <- function(data, flank=2) 
# This function computes the count sum for the bin pairs that are neighbouring
# each specified bin pair in 'data'. For inter-chromosomal bin pairs, it does
# so by considering a cross-like shape, centered on the current bin pair.  For
# intra-chromosomal bin pairs, it does so by considering adjacent bin pairs on
# the same diagonal. Note that for this to work, you need 'filter=1L' during
# count loading. It's a bit of a kludge but it's the simplest way to do it.
#
# written by Aaron Lun
# Created 23 April 2014
# Modified 3 November 2014
{
	flank <- as.integer(flank)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	rdata <- .delimitFragments(regions(data))
	last.id <- rdata$end
	first.id <- rdata$start
	names(last.id) <- names(first.id) <- rdata$chr

	np <- nrow(data)
	by.chr <- split(1:np, as.character(seqnames(anchors(data))))
	output <- matrix(0L, ncol=ncol(data), nrow=np)
   	out.n <- integer(np)

	# Running through each pair of chromosomes.
	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, as.character(seqnames(targets(data[next.chr,]))))
		a.len <- last.id[[anchor]] - first.id[[anchor]] + 1L

		for (target in names(next.chr)) {
			current.pair <- next.chr[[target]]
			all.a <- data@anchor.id[current.pair] - first.id[[anchor]] 
			all.t <- data@target.id[current.pair] - first.id[[target]]
			rel.counts <- data@counts[current.pair,,drop=FALSE]

			if (target==anchor) {
				# Doing diagonal-based background estimation.
				diagonal <- all.a - all.t
				position <- all.t
				o <- order(diagonal, position)
				collected <-.Call(cxx_collect_diagonal, diagonal[o], position[o], rel.counts[o,], flank, a.len)
				if (is.character(collected)) { stop(collected) }

				collected[[1]][o,] <- collected[[1]]
				collected[[2]][o] <- collected[[2]]
				all.counts <- collected[[1]] - rel.counts
				all.n <- collected[[2]] - 1L
			} else {
				t.len <- last.id[[target]] - first.id[[target]] + 1L

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
