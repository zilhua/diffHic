compartmentalize <- function(data, centers=2, inter=FALSE, ...)
# Computes compartments for every chromosome, using the intra-chromosomal
# contact maps that have been corrected for distance effects.
#
# written by Aaron Lun
# created 26 May 2015
{
	if (!inter) {
		is.intra <- !is.na(getDistance(data))
		trended <- filterTrended(data[is.intra,])
		resids <- rep(NA, length(is.intra))
		resids[is.intra] <- trended$abundance - trended$threshold
	} else { 
		trended <- filterTrended(data)
		resids <- trended$abundance - trended$threshold
	}
	dist2trend <- approxfun(x=trended$log.distance, y=trended$threshold, rule=2)

	if (!inter) { 
		stored <- list()
		for (chr in seqlevels(regions(data))) { 
			mat <- as.matrix(data, first=chr, fill=resids)
			if (any(dim(mat)==0L)) { next }
			mat <- .fillZeros(mat, data, dist2trend)
			mat <- .debiasBins(mat)
					
			# Actually partitioning by location.
			comp <- .safeCluster(mat, centers=centers, ...)
			names(comp) <- rownames(mat)
			stored[[chr]] <- list(compartment=comp, matrix=mat)
		}
	} else {
		# Using the entirety of the interaction space.				
		mat <- as.matrix(data, fill=resids)
		mat <- .fillZeros(mat, data, dist2trend)
		mat <- .debiasBins(mat)
		comp <- .safeCluster(mat, centers=centers, ...)
		stored <- list(compartment=comp, matrix=mat)
	}

	return(stored)
}

.fillZeros <- function(mat, data, dist2trend) {
	# Filling NA's (i.e., zero's). Using mid-distance, interpolating to get the trend.
	lost <- which(is.na(mat), arr.ind=TRUE)
	lost.dist <- abs(mid(ranges(anchors(data[lost[,1]]))) - mid(ranges(targets(data[lost[,2]]))))
	lost.dist <- log10(lost.dist + exptData(data)$width)
	mat[is.na(mat)] <- .makeEmpty(data) - dist2trend(lost.dist)
	return(mat)
}

.safeCluster <- function(mat, centers, ...) { 
	if (nrow(mat) > centers) { 
		out <- kmeans(mat, centers=centers, ...)
		comp <- out$cluster
	} else {
		comp <- 1:nrow(mat)
	}
	return(comp)
}

.debiasBins <- function(mat) {
	# Just subtracting half from both rows and columns. This is equivalent to 
	# dividing by square root of coverage, which works pretty well in place of
	# a more rigorous iterative approach (check out Rao's supplementaries).
	rwm <- log2(rowMeans(2^mat))
	mat <- mat - rwm/2
	mat <- t(t(mat) - rwm/2)
	return(mat)
}

