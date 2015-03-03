clusterPairs <- function(data, tol, upper=1e6) 
# This function examines the bin pairs in two-dimensional space and 
# clusters those pairs which are close together. Specifically, it does 
# so if the Chebyshev distance between the regions is less than 'tol'.
# It also splits clusters which are greater than 'upper' by partitioning
# them into smaller clusters of equal size.
#
# written by Aaron Lun
# created 6 December 2013
# last modified 3 March 2015
{
	region <- regions(data)
	allchrs <- as.character(seqnames(region))
	achrs <- allchrs[data@anchors]
	tchrs <- allchrs[data@targets]
	tol <- as.integer(tol)
	stopifnot(tol>=0L) # Minimum overlap not supported.
	upper <- as.integer(upper)

	# Figuring out which are blocks of separate chromosomes.	
	ro <- order(achrs, tchrs)
	achrs <- achrs[ro]
	tchrs <- tchrs[ro]
	n <- length(ro)
	is.new <- which(c(TRUE, achrs[-1]!=achrs[-n] | tchrs[-1]!=tchrs[-n]))
	upnext <- c(is.new[-1]-1L, n)
	
	# Setting up the starts and ends.
	astarts <- start(region[data@anchors[ro]])
	aends <- end(region[data@anchors[ro]])+1L
	tstarts <- start(region[data@targets[ro]])
	tends <- end(region[data@targets[ro]])+1L

	# Now, running through.
	all.ids <- integer(n)
	bonus <- 0L
	for (i in 1:length(upnext)) {
		current <- is.new[i]:upnext[i]			
		curas <- astarts[current]	
		curae <- aends[current]	
		curts <- tstarts[current]	
		curte <- tends[current]
		po <- order(curas, curts)
	
		out <- .Call(cxx_cluster_2d, curas[po], curts[po], curae[po], curte[po], tol)
		if (is.character(out)) { stop(out) }
		out[po] <- out
		if (length(upper)) { 
			out <- .Call(cxx_split_clusters, out, curas, curts, curae, curte, upper)
			if (is.character(out)) { stop(out) }
			xo <- order(out)
			out[xo]<-cumsum(c(TRUE, diff(out[xo])!=0L)) # Cleaning it up a little, to keep the IDs reasonably tight.
		}

		all.ids[current] <- out + bonus
		bonus <- bonus + max(out)
	}

	all.ids[ro] <- all.ids
	return(all.ids)
}
