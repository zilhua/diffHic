filterDirect <- function(data, ...)
# Implements the direct filtering method on the abundances of 
# inter-chromosomal bin pairs.
#
# written by Aaron Lun
# created 5 March 2015
{
	all.chrs <- seqnames(regions(data))
	is.inter <- as.logical(all.chrs[anchors(data, id=TRUE)]!=all.chrs[targets(data, id=TRUE)])
	ave.ab <- scaledAverage(asDGEList(data[is.inter,]), ...)
	threshold <- .getInterThreshold(all.chrs, ave.ab)
	return(list(abundances=ave.ab, threshold=threshold))
}

.getInterThreshold <- function(all.chrs, inter.ab) { 
	# Getting the total number of inter-chromosomal bins.
	n.bins <- as.numeric(runLength(all.chrs))
	total.bins <- sum(n.bins)
	n.inter <- total.bins * (total.bins + 1L)/2L - sum(n.bins * (n.bins + 1L)/2L)
	prop.kept <- length(inter.ab)/n.inter

	# Getting the threshold.
	if (prop.kept >= 1) { 
		threshold <- median(inter.ab) 
	} else if (prop.kept < 0.5) { 
		threshold <- scaledAverage(DGEList(0, lib.size=mean(data$totals)), ...) 
	} else { 
		threshold <- quantile(inter.ab, 1-0.5/prop.kept)
		names(threshold) <- NULL
	}
	return(threshold)
}

filterTrended <- function(data, span=0.25, ...)
# Implements the trended filtering method on the abundances of 
# inter-chromosomal bin pairs. Don't bother getting the missing
# entries, they'll be ignored due to robustification in loess.
#
# written by Aaron Lun
# created 5 March 2015
{
	dist <- getDistance(data, type="mid")
	log.dist <- log10(dist + exptData(data)$width)
	ave.ab <- scaledAverage(asDGEList(data), ...)
	trend.threshold <- loessFit(x=log.dist, y=ave.ab, span=span)$fitted

	# Using the direct threshold.
	direct.threshold <- .getInterThreshold(seqnames(regions(data)), ave.ab[is.na(log.dist)])
	trend.threshold[is.na(log.dist)] <- direct.threshold
	return(list(abundances=ave.ab, threshold=trend.threshold, log.distance=log.dist)) 
}

