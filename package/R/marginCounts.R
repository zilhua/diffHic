marginCounts <- function(files, fragments, width=500000, restrict=NULL)
# Gets the marginal counts i.e. sum of counts for each fragment.
{
	width <- as.integer(width)
    if (width < 0) { stop("width must be a non-negative integer") }
    new.pts <- .getBinID(fragments, width)
	chrs <- seqlevels(fragments)
	total.bins <- length(new.pts$region)
	stopifnot(max(new.pts$id)==total.bins) 

    # Running through each pair of chromosomes.
	nlibs <- length(files)
    total.counts <- matrix(0L, length(new.pts$region), nlibs)
	full.sizes <- integer(nlibs)

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
	keep <- rowSums(total.counts)>0.5
	return(list(counts=total.counts[keep,,drop=FALSE], totals=full.sizes, region=new.pts$region[keep]))
}

