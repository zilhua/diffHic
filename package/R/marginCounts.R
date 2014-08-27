marginCounts <- function(files, fragments, width=500000, restrict=NULL)
# Gets the marginal counts i.e. sum of counts for each bin or region.
# This is useful to determine the `genomic coverage' of each region,
# based on the number of Hi-C read pairs involving that region.
#
# written by Aaron Lun
# Some time ago.
{
	width <- as.integer(width)
    if (width < 0) { stop("width must be a non-negative integer") }
    new.pts <- .getBinID(fragments, width)
	total.bins <- length(new.pts$region)
	stopifnot(max(new.pts$id)==total.bins) 

	nlibs <- length(files)
    total.counts <- matrix(0L, length(new.pts$region), nlibs)
	full.sizes <- integer(nlibs)
	chrs <- seqlevels(fragments)

    # Running through each pair of chromosomes.
    overall <- .loadIndices(files)
    for (anchor in names(overall)) {
        stopifnot(anchor %in% chrs)
		if (!is.null(restrict) && !(anchor %in% restrict)) { next }
		current <- overall[[anchor]]
		for (target in names(current)) {
			stopifnot(target %in% chrs)
			if (!is.null(restrict) && !(target %in% restrict)) { next }
    
      		pairs <- .baseHiCParser(current[[target]], files, anchor, target)
           	full.sizes <- full.sizes+sapply(1:length(pairs), FUN=function(x) { sum(pairs[[x]]$count) })

			# Aggregating them in C++ to get the count combinations and location of each bin.
            out <- .Call(cxx_count_marginals, pairs, new.pts$id, total.bins)
            if (is.character(out)) { stop(out) }
			total.counts <- total.counts+out
		}
	}
	
	# Aggregating all elements.
	retained <- which(rowSums(total.counts)>0.5)
	return(.DIList(counts=total.counts[retained,,drop=FALSE], totals=full.sizes, 
			anchors=retained, targets=retained, regions=new.pts$region))
}

