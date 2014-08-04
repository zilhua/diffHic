getDistance <- function(data, type=c("mid", "gap", "span")) 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
#
# written by Aaron Lun
# 22 April, 2014
{
	is.same <- as.logical(seqnames(data$region[data$pairs[,1]])==seqnames(data$region[data$pairs[,2]]))
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

getArea <- function(data)
# Computing the number of restriction fragment pairs in the interaction space.
# This allows adjustment of abundances for comparison between differently-sized areas.
# Coercion to double is necessary in case the number of fragments overflows the integer type.
# 
# written by Aaron Lun
# 30 July, 2014
{
#	is.same <- data$pairs[,1] == data$pairs[,2]
	curnfrag <- as.double(data$region$nfrags[data$pairs[,1]])
	returned <- curnfrag * data$region$nfrags[data$pairs[,2]]
#	returned[is.same] <- curnfrag[is.same]*(curnfrag[is.same]+1)/2
	returned
}

# Special behaviour is necessary on the diagonal, as only half the fragments are actually used.
# Note that this won't work on partially overlapping ranges, because there's not enough
# information to determine how many restriction fragments are overlapped. So, I just didn't
# put it in at all for simplicity.

