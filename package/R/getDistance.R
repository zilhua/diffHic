getDistance <- function(data, type=c("mid", "gap", "span")) 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
#
# written by Aaron Lun
# created 22 April 2014
# last modified 3 March 2015
{
	all.as <- anchors(data)
	all.ts <- targets(data)
	is.same <- as.logical(seqnames(all.as)==seqnames(all.ts))
	all.as <- all.as[is.same]
	all.ts <- all.ts[is.same]
	output <- rep(NA, nrow(data))

	type <- match.arg(type)
	if (type=="gap") {
		output[is.same] <- pmax(start(all.as), start(all.ts)) - pmin(end(all.as), end(all.ts)) + 1
	} else if (type=="span") {
		output[is.same] <- pmax(end(all.as), end(all.ts)) - pmin(start(all.as), start(all.ts)) + 1
	} else if (type=="mid") {
		output[is.same] <- mid(ranges(all.as)) - mid(ranges(all.ts))
	}
	return(output)
}

getArea <- function(data, bp=TRUE)
# Computes the area of the interaction space, either in base pair-squared terms
# or in terms of pairs of restriction fragments. This allows adjustment of
# abundances for comparison between differently-sized areas. Special behaviour
# is necessary on the diagonal, as reflection halves the space. Coercion to 
# double is necessary to prevent overflows of the integer type.
# 
# written by Aaron Lun
# created 30 July 2014
# last modified 5 March 2015
{
	ax <- data@anchors
	tx <- data@targets

	if (bp) {
		cur.width <- as.double(width(data@regions))
		returned <- cur.width[ax] * cur.width[tx]

		# Accounting for special behaviour around the diagonal.	It you don't halve,
		# you'll have to double every other (unreflected) area.
		overlap <- getDistance(data, type="gap")
		is.olap <- !is.na(overlap) & overlap < -0.5
		lap.dist <- abs(overlap)[is.olap]
		self.lap.area <- lap.dist * (lap.dist - 1)/2		
		returned[is.olap] <- returned[is.olap] - self.lap.area
	} else {
		is.same <- ax==tx
		curnfrag <- as.double(data@regions$nfrags[ax])
		returned <- curnfrag * data@regions$nfrags[tx]
		returned[is.same] <- curnfrag[is.same]*(curnfrag[is.same]+1)/2

		# Detour to protect against overlapping regions.
		left.edge <- pmax(start(data@regions)[ax], start(data@regions)[tx])
		right.edge <- pmin(end(data@regions)[ax], end(data@regions)[tx])
		is.partial <- !is.same & right.edge >= left.edge & 
			as.logical(seqnames(data@regions)[ax]==seqnames(data@regions)[tx]) 

		if (any(is.partial)) { 
			right.edge <- right.edge[is.partial]
			left.edge <- left.edge[is.partial]
			by.chr <- split(1:sum(is.partial), as.character(seqnames(data@regions)[ax][is.partial]))
			fragments <- exptData(data)$param$fragments
			fdata <- .delimitFragments(fragments)

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
	}

	return(returned)
}
