totalCounts <- function(files, fragments=NULL, restrict=NULL)
# This function gets the total counts in a bunch of files.  This is designed
# for whenever the total counts must be rapidly extracted, without the need to
# count across the interaction space.
#
# written by Aaron Lun
# 17 September 2014
{
	nlibs <- length(files)
	if (nlibs==0) { 
		stop("number of libraries must be positive")
	} 
	if (!is.null(fragments)) { chrs <- seqlevels(fragments) }
	full.sizes <- integer(nlibs)

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files)
    for (anchor in names(overall)) {
		if (!is.null(fragments)) { stopifnot(anchor %in% chrs) }
	    if (!is.null(restrict) && !(anchor %in% restrict)) { next }
        current <- overall[[anchor]]

		for (target in names(current)) {
			if (!is.null(fragments)) { stopifnot(target %in% chrs) }
	        if (!is.null(restrict) && !(target %in% restrict)) { next }

			# Getting totals.
			pairs <- .baseHiCParser(current[[target]], files, anchor, target)
			for (lib in 1:length(pairs)) { full.sizes[lib] <- full.sizes[lib] + sum(pairs[[lib]]$count) }
		}
	}

	return(full.sizes)
}


