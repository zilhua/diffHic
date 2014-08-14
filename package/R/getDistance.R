getDistance <- function(data, type=c("mid", "gap", "span")) 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
#
# written by Aaron Lun
# 22 April, 2014
{
	is.same <- as.logical(seqnames(data$region)[data$pairs[,1]]==seqnames(data$region)[data$pairs[,2]])
	all.as <- data$region[data$pairs[is.same,1]]
	all.ts <- data$region[data$pairs[is.same,2]]
	output <- rep(NA, nrow(data$pairs))

	type <- match.arg(type)
	if (type=="gap") {
		output[is.same] <- pmin(end(all.as), end(all.ts)) - pmax(start(all.as), start(all.ts))
	} else if (type=="span") {
		output[is.same] <- pmax(end(all.as), end(all.ts)) - pmin(start(all.as), start(all.ts))
	} else if (type=="mid") {
		output[is.same] <- mid(ranges(all.as)) - mid(ranges(all.ts))
	}
	return(output)
}

getarea <- function(data, fragments=null)
# Computing the number of restriction fragment pairs in the interaction space.
# This allows adjustment of abundances for comparison between differently-sized areas.
# Special behaviour is necessary on the diagonal, as only half the fragments are actually used.
# 'fragments' should be specified if there is any risk of partial overlaps.
# Coercion to double is necessary in case the number of fragments overflows the integer type.
# 
# written by Aaron Lun
# 30 July, 2014
# modified 14 August 2014
{
	ax <- data$pairs[,1]
	tx <- data$pairs[,2]	
	is.same <- ax==tx
	curnfrag <- as.double(data$region$nfrags[ax])
	returned <- curnfrag * data$region$nfrags[tx]
	returned[is.same] <- curnfrag[is.same]*(curnfrag[is.same]+1)/2

	if (!is.null(fragments)) { 
		# Detour to protect against overlapping regions.
		fdata <- .checkFragments(fragments)

		left.edge <- pmax(start(data$region)[ax], start(data$region)[tx])
		right.edge <- pmin(end(data$region)[ax], end(data$region)[tx])
		is.partial <- !is.same & right.edge >= left.edge & 
			as.logical(seqnames(data$region)[ax]==seqnames(data$region)[tx]) 

		right.edge <- right.edge[is.partial]
		left.edge <- left.edge[is.partial]
		by.chr <- split(1:sum(is.partial), as.character(seqnames(data$region)[ax][is.partial]))

		for (x in 1:length(fdata$chr)) {
			current.chr <- fdata$chr[x]
			curdex <- by.chr[[current.chr]]
			if (is.null(curdex)) { next }
		
			indices <- fdata$start[x]:fdata$end[x]
			right.olap <- match(right.edge[curdex], end(fragments)[indices])
			left.olap <- match(left.edge[curdex], start(fragments)[indices])
 	    	if (any(is.na(right.olap)) || any(is.na(left.olap))) { stop("region boundaries should correspond to restriction fragment boundaries") }
		
			n.overlap <- right.olap - left.olap + 1	
			returned[is.partial][curdex] <- returned[is.partial][curdex] - n.overlap*(n.overlap-1)/2
		}
	}

	return(returned)
}
