countNeighbors <- function(data, flank, type=c("both", "anchor", "target")) 
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
	rdata <- .checkFragments(data$region)
	last.id <- rdata$end
	first.id <- rdata$start
	names(last.id) <- names(first.id) <- rdata$chr

	stopifnot(all(data$pairs$anchor.id >= data$pairs$target.id))
	rownames(data$pairs) <- NULL
	by.chr <- split(data$pairs, as.character(seqnames(data$region[data$pairs$anchor.id])))
	np <- nrow(data$pairs)
	output <- matrix(0L, ncol=ncol(data$counts), nrow=np)
   	out.n <- numeric(np)

	# Running through each pair of chromosomes.
	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, as.character(seqnames(data$region[next.chr$target.id])))
		if (!anchor %in% names(last.id)) { stop("chromosome missing in the fragment list") }

		for (target in names(next.chr)) {
			cur.chr <- next.chr[[target]]

			# We need to double it, to account for reflection around the diagonal. 
			# We obviously don't need to add the diagonal points twice.
			on.diag <- cur.chr$anchor.id==cur.chr$target.id
			reflected <- !on.diag & (cur.chr$target.id + flank >= cur.chr$anchor.id)
			combined.anchor <- c(cur.chr$anchor.id, cur.chr$target.id[reflected])
			combined.target <- c(cur.chr$target.id, cur.chr$anchor.id[reflected])
			keep <- 1:nrow(cur.chr)
			relevant <- as.integer(rownames(cur.chr))
			rel.counts <- data$counts[relevant,,drop=FALSE]
			combined.counts <- rbind(rel.counts, rel.counts[reflected,,drop=FALSE])

			# Computing the local background, for fixed anchor or target.
			all.counts <- all.n <- 0L
			if (type!="anchor") { 
				o <- order(combined.anchor, combined.target)
				collected <- .Call(cxx_collect_background, combined.anchor[o], combined.target[o], combined.counts[o,], flank)
				if (is.character(collected)) { stop(collected) }
				collected[o,] <- collected
				all.counts <- collected[keep,] - rel.counts

				upper <- pmin(cur.chr$target.id+flank, last.id[[target]]) 
				lower <- pmax(cur.chr$target.id-flank, first.id[[target]]) 
				all.n <- upper - lower
			}
			if (type!="target") { 
				o <- order(combined.target, combined.anchor)
				collected <- .Call(cxx_collect_background, combined.target[o], combined.anchor[o], combined.counts[o,], flank)
				if (is.character(collected)) { stop(collected) }
				collected[o,] <- collected
				all.counts <- all.counts + collected[keep,] - rel.counts
				
				upper <- pmin(cur.chr$anchor.id+flank, last.id[[anchor]]) 
				lower <- pmax(cur.chr$anchor.id-flank, first.id[[anchor]]) 
				all.n <- all.n + upper - lower
			}

			# Points right on the diagonal count the neighborhood twice upon reflection, so division is necessary.
			if (type=="both") { 
				all.counts[on.diag,] <- all.counts[on.diag,]/2L
				all.n[on.diag] <- all.n[on.diag]/2L
			}
			output[relevant,] <- all.counts
			out.n[relevant] <- all.n
		}
	}

	# Returning the collected counts.
	return(list(counts=output, n=out.n))
}
