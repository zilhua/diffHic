filterPeaks <- function(data, flank=5, prior.count=2, abundances=NULL)
# This function identifies the highest-abundance neighbour in the interaction space
# for each bin pair in `data`. The aim is to compare the abundance of each element
# with the abundance of its neighbour. 
#
# written by Aaron Lun
# created 23 April 2014
# last modified 23 March 2015
{
	flank <- as.integer(flank)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	rdata <- .delimitFragments(regions(data))
	last.id <- rdata$end
	first.id <- rdata$start
	names(last.id) <- names(first.id) <- rdata$chr

	all.ab <- abundances
	if (is.null(all.ab)) { all.ab <- aveLogCPM(asDGEList(data), prior.count=0) }

	# Rescaling to count-level data with at least 6 dp, for stable calculations with integers.
	# This should be enough precision while avoiding overrun of the integer type.
	scaling <- log2(mean(data$totals)/1e6 * ncol(data))
	MULT <- 1e6
	back2count <- function(ab) { 
		newval <- 2^(ab + scaling)
		lower <- as.integer(newval)
		list(val=newval, int=lower, dec=as.integer(round((newval-lower)*MULT)))
	}
	prior.count <- prior.count * ncol(data)	
	
	# Running through each pair of chromosomes.
	np <- nrow(data)
	by.chr <- split(1:np, as.character(seqnames(anchors(data))))
	aid <- anchors(data, id=TRUE)
	tid <- targets(data, id=TRUE)
	output <- numeric(np)

	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, as.character(seqnames(targets(data[next.chr,]))))
		a.len <- last.id[[anchor]] - first.id[[anchor]] + 1L

		for (target in names(next.chr)) {
			current.pair <- next.chr[[target]]
			all.a <- aid[current.pair] - first.id[[anchor]] 
			all.t <- tid[current.pair] - first.id[[target]]
			t.len <- last.id[[target]] - first.id[[target]] + 1L

			rel.ab <- all.ab[current.pair]
			converted <- back2count(rel.ab) 

			# Using the quadrant with the maximum average.
			o <- order(all.a, all.t)
			collected <- .Call(cxx_quadrant_bg, all.a[o], all.t[o], 
				converted$int[o], converted$dec[o], MULT, 
				flank, a.len, t.len, anchor==target)
			if (is.character(collected)) { stop(collected) }
			collected[o] <- collected

			output[current.pair] <- log2((converted$val+prior.count)/(collected+prior.count))
		}
	}

	# Returning the collected log-fold-changes.
	return(output)
}

# Note that, when trend!="none" and we're looking on the same chromosome,
# rescaled distance-adjusted values won't be on the same scale as the original
# counts. This isn't a major problem as we're calculating the relative
# enrichment anyway. However, there mightn't be enough decimal places to cover
# distance-adjusted values that are very low. Hopefully, more values should be
# around 1, so it shouldn't be a problem.

