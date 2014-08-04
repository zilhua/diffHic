countArea <- function(files, fragments, region1, region2, self.circles=FALSE, type="any") 
# This counts the number of read pairs mapped between each of 'region1' and 'region2'.
# This is very similar to connectCounts, except that the C++ counting function is 
# silghtly different.
#
# written by Aaron Lun
# 25 July 2014
{
	fdata <- .checkFragments(fragments)

	# Finding the overlaps. This tries to reduce the amount of work by only
	# running findOverlaps on unique 'region1' and 'region2'.
	nr <- length(region1)
	stopifnot(nr!=length(region2))
	olap.1 <- .retrieveHits(region1, fragments, type=type)
	renewed.1 <- olap.1$hits
	start.1 <- olap.1$start
	end.1 <- olap.1$end
	olap.2 <- .retrieveHits(region2, fragments, type=type)
	renewed.2 <- olap.2$hits
	start.2 <- olap.2$start
	end.2 <- olap.2$end

#	o.1 <- order(region1)
#	o.2 <- order(region2)
#	oregion1 <- region1[o.1]
#	oregion2 <- region2[o.2]
#
#	u.1 <- !duplicated(oregion1)
#	u.2 <- !duplicated(oregion2)
#	olap.1 <- findOverlaps(fragments, oregion1[u.1], type=type)
#	olap.2 <- findOverlaps(fragments, oregion2[u.2], type=type)
#
#	ulen.1 <- cumsum(u.1)
#	uu1 <- split(o.1, ulen.1)
#	ulen.1 <- rle(ulen.1)$length
#	ulen.2 <- cumsum(u.2)
#	uu2 <- split(o.2, ulen.2)
#	ulen.2 <- rle(ulen.2)$length
#	
#	renewed.1 <- unlist(uu1[subjectHits(olap.1)])
#	renewed.2 <- unlist(uu2[subjectHits(olap.2)])
#	rle.1 <- rle(rep(queryHits(olap.1), ulen.1))
#	end.1 <- cumsum(rle.1$length) + 1L
#	start.1 <- end.1 - rle.1$length 
#	rle.2 <- rle(rep(queryHits(olap.2), ulen.2))
#	end.2 <- cumsum(rle.2$length) + 1L
#	start.2 <- end.2 - rle.2$length 

	# Picking out the combinations, so we don't have to load everything.
	ac <- as.character(seqnames(region1))
	tc <- as.character(seqnames(region2))
	a.is.anchor <- match(ac, fdata$chrs) >= match(tc, fdata$chrs)
	all.ac <- ifelse(a.is.anchor, ac, tc)
	all.tc <- ifelse(a.is.anchor, tc, ac)
	o <- order(all.ac, all.tc)
	all.ac <- all.ac[o]
	all.tc <- all.tc[o]
	is.diff <- c(1L, which(all.ac[-1]!=all.ac[-nr] | all.tc[-1]!=all.tc[-nr])+1L)

	# Pulling out the reads and finding out what they line up with.
	nlib <- length(files)
	collected <- matrix(0L, nr, nlib)
	for (i in is.diff) {
		current.a <- all.ac[i]
		current.t <- all.tc[i]
		for (j in 1:nlib) {
			current <- .getPairs(files[j], current.a, current.t) 
			out <- .Call(cxx_count_area, current,
				start.1, end.1, renewed.1,
				start.2, end.2, renewed.2)
			if (is.character(out)) { stop(out) }
			collected[,j] <- collected[,j] + out
		}
	}
	
	return(collected)
}
