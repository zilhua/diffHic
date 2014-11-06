enrichedGap <- function(data, bin.size, mode=15, flank=2, min.flank=0, 
	per.chr=TRUE, prior.count=5, span=0.3, ...) 
# This function identifies the highest-abundance neighbour in the interaction space
# for each bin pair in `data`. The aim is to compare the abundance of each element
# with the abundance of its neighbour. Some adjustment is necessary to account for
# the effect of distance on abundance when comparing between distances.
#
# written by Aaron Lun
# Created 23 April 2014
# Modified 3 November 2014
{
	flank <- as.integer(flank)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	min.flank <- as.integer(min.flank)
	if (min.flank < 0L || min.flank >= flank) { stop("minimum flank must be non-negative and less than the flank width") }
	rdata <- .delimitFragments(regions(data))
	last.id <- rdata$end
	first.id <- rdata$start
	names(last.id) <- names(first.id) <- rdata$chr

	# Checking the mode.
	mode <- as.integer(mode)
	if (mode <= 0L || mode >= 16L) { stop("mode must be a positive integer under 16") }
	vertical <- bitwAnd(mode, 0x1) != 0L
	horizontal <- bitwAnd(mode, 0x2) != 0L
	diag.1 <- bitwAnd(mode, 0x4) != 0L
	diag.2 <- bitwAnd(mode, 0x8) != 0L

	# Instantiating variables.
	all.ab <- aveLogCPM(asDGEList(data), prior.count=prior.count)
	if (missing(bin.size)) { 
		warning("guessing bin size from median region width")
		bin.size <- median(width(regions(data)))
	}
	log.dist <- log10(getDistance(data) + bin.size)
	if (!per.chr) { fit <- loessFit(y=all.ab, x=log.dist, span=span) }

	np <- nrow(data)
	by.chr <- split(1:np, as.character(seqnames(anchors(data))))
	output <- numeric(np)
	null.ab <- aveLogCPM(rbind(integer(ncol(data))), lib.size=totals(data), prior.count=prior.count)

	# Running through each pair of chromosomes.
	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, as.character(seqnames(targets(data[next.chr,]))))
		a.len <- last.id[[anchor]] - first.id[[anchor]] + 1L

		for (target in names(next.chr)) {
			current.pair <- next.chr[[target]]
			all.a <- data@anchor.id[current.pair] - first.id[[anchor]] 
			all.t <- data@target.id[current.pair] - first.id[[target]]
			rel.ab <- all.ab[current.pair]
			t.len <- last.id[[target]] - first.id[[target]] + 1L

			if (target==anchor) {
				# Adjusting abundances for distance. 
				if (per.chr) {
					if (length(current.pair)>1L) { 
						adj <- loessFit(y=rel.ab, x=log.dist[current.pair], span=span)$fitted
					} else {
						adj <- 0 
					}
				} else {
					adj <- fit$fitted[current.pair]
				}
				rel.ab <- rel.ab -adj
			} 

			# Computing neighbourhood for same anchor, different target.
			max.ab <- -Inf
			if (horizontal) { 
				o <- order(all.a, all.t)
				collected <- .Call(cxx_max_background, all.a[o], all.t[o], rel.ab[o], flank, min.flank)
				if (is.character(collected)) { stop(collected) }
				collected[o] <- collected
				max.ab <- pmax(max.ab, collected)
			}
			if (vertical) { 
				o <- order(all.t, all.a)
				collected <- .Call(cxx_max_background, all.t[o], all.a[o], rel.ab[o], flank, min.flank)
				if (is.character(collected)) { stop(collected) }
				collected[o] <- collected
				max.ab <- pmax(max.ab, collected)
			}
			if (diag.1) {
				d.up <- all.a - all.t
				o <- order(d.up, all.t)
				collected <- .Call(cxx_max_background, d.up[o], all.t[o], rel.ab[o], flank, min.flank)
				if (is.character(collected)) { stop(collected) }
				collected[o] <- collected
				max.ab <- pmax(max.ab, collected)
			}
			if (diag.2) {
				d.down <- all.a + all.t
				o <- order(d.down, all.t)
				collected <- .Call(cxx_max_background, d.down[o], all.t[o], rel.ab[o], flank, min.flank)
				if (is.character(collected)) { stop(collected) }
				collected[o] <- collected
				max.ab <- pmax(max.ab, collected)
			}
			
			max.ab[is.infinite(max.ab)] <- null.ab
			output[current.pair] <- rel.ab - max.ab
		}
	}

	# Returning the collected counts.
	return(output)
}
