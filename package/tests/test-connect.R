######################################################################################
# This tests the functionality of connectCounts.

suppressPackageStartupMessages(require(diffHic))
suppressPackageStartupMessages(require(rhdf5))
chromos <- c(chrA=100, chrB=80)
source("simcounts.R")

simranges <- function(cuts, nranges, min.size=1000, max.size=10000)
# Generates simulated ranges.
{
    ranges <- list()
	for (chr in seqlevels(cuts)) {
		chr.len <- seqlengths(cuts)[[chr]] 
		max.val <- chr.len - min.size
		range.start <- round(runif(nranges, 1, max.val))
		range.end <- pmin(chr.len, round(range.start + runif(nranges, min.size, max.size)))
		ranges[[chr]] <- GRanges(chr, IRanges(range.start, range.end))
	}		
	names(ranges) <- NULL
	suppressWarnings(ranges <- do.call(c, ranges))
	return(ranges)	
}

reconstruct <- function(pairs, counts) {
	if (is.null(pairs) || nrow(pairs)<=1L) { return(list(pairs=pairs, counts=counts)) }
	counts <- as.matrix(counts)
	o <- order(pairs[,1], pairs[,2])
	pairs <- pairs[o,]
	for (i in 1:ncol(counts)) { counts[,i] <- cumsum(counts[o,i]) }
	last.diff <- c(diff(pairs[,1])!=0L | diff(pairs[,2])!=0L, TRUE)
	my.count <- apply(counts, 2, FUN=function(x) { diff(c(0L, x[last.diff])) })
	if (is.null(dim(my.count))) { my.count <- rbind(my.count) } 
	return(list(pairs=pairs[last.diff,,drop=FALSE], counts=my.count))
}

refline <- function(dirs, cuts, ranges, filter=20L, type="any") {
	everypair <- everycount <- list()
	o <- order(ranges)
	ranges <- ranges[o]

	# Determining which ranges each restriction fragment overlaps.
	cur.olap <- findOverlaps(cuts, ranges, type=type)
	cur.rle <- rle(queryHits(cur.olap))
	cur.end <- cumsum(cur.rle$length)
	cur.start <- cur.end - cur.rle$length + 1L
	cur.hits <- subjectHits(cur.olap)

	for (d in 1:length(dirs)) {
		allpairs <- allcounts <- list()
    	x <- h5ls(dirs[d])
		x <- x[x$otype=="H5I_DATASET",]

	    for (k in 1:length(chromos)) {
	        cur.k<-names(chromos)[k]
	        for (l in 1:k) {
	            cur.l<-names(chromos)[l]
				if (!any(basename(x$group)==cur.k & x$name==cur.l)) { next }
				counts <- h5read(dirs[d], file.path(cur.k, cur.l))
				for (xx in 1:ncol(counts)) { attributes(counts[,xx]) <- NULL }

				# Need in both.
				collected <- list()
				collected.counts <- list()
				matched.a <- match(counts$anchor.id, cur.rle$values)
				matched.t <- match(counts$target.id, cur.rle$values)
				in.both <- !is.na(matched.a) & !is.na(matched.t)

				# Determining which ranges each pair overlaps.
				for (j in which(in.both)) {
					ja <- matched.a[j]
					jt <- matched.t[j]
					in.a <- cur.hits[cur.start[ja]:cur.end[ja]]
					in.t <- cur.hits[cur.start[jt]:cur.end[jt]]
					additionals <- as.matrix(expand.grid(in.a, in.t))
					flipped <- additionals[,2] >= additionals[,1]
					additionals[flipped,] <- additionals[flipped,2:1]

					additionals <- reconstruct(additionals, 1:nrow(additionals))$pairs # Eliminating redundant elements for each pair.
					idex <- length(collected) + 1L
					collected[[idex]] <- additionals
					collected.counts[[idex]] <- rep(counts$count[j], nrow(additionals))
				}
			
				# Assembling summary counts for this chromosome combination in this library.
				collected <- do.call(rbind, collected)
				collected.counts <- do.call(c, collected.counts)
				if (is.null(collected.counts)) { next }
				out <- reconstruct(collected, collected.counts)
				idex <- length(allpairs) + 1L
				allpairs[[idex]] <- out$pairs
				allcounts[[idex]] <- out$counts
			}
		}

		# No need to summarize here, combinations will be different between chromosome pairs.
		allpairs <- do.call(rbind, allpairs)
		allcounts <- unlist(allcounts)
		actually <- matrix(0L, ncol=length(dirs), nrow=length(allcounts))
		actually[,d] <- allcounts
		idex <- length(everypair)
		everypair[[idex+1L]] <- allpairs
		everycount[[idex+1L]] <- actually
	}

	# Redefining regions.
	so <- subjectHits(cur.olap)
	qo <- queryHits(cur.olap)
	new.starts <- by(start(cuts[qo]), INDICES=so, FUN=min)
	new.ends <- by(end(cuts[qo]), INDICES=so, FUN=max)
	new.num <- by(start(cuts[qo]), INDICES=so, FUN=length)
	acquired <- as.integer(names(new.starts))
	
	ranges2 <- ranges
	start(ranges2)[acquired] <- as.integer(new.starts)
	end(ranges2)[acquired] <- as.integer(new.ends)
	full.num <- integer(length(ranges2))
	full.num[acquired] <- as.integer(new.num)
	ranges2$nfrags <- full.num
	ranges2$original <- o

	# Aggregating results between libraries.
	everypair <- do.call(rbind, everypair)
	everycount <- do.call(rbind, everycount)
	if (is.null(everycount) || nrow(everycount)==0L) { 
		final <- list(pairs=data.frame(anchor.id=integer(0), target.id=integer(0)), 
				counts=matrix(0L, ncol=length(dirs), nrow=0), region=ranges2)
		return(final)
	}
	final <- reconstruct(everypair, everycount)
	keep <- rowSums(final$counts) >= filter

	# Determining which one is anchor or target.
	left <- final$pairs[keep,1]
	right <- final$pairs[keep,2]
	matched <- match(as.character(seqnames(ranges)), runValue(seqnames(cuts)))
	rank <- integer(length(ranges))
	rank[order(matched, start(ranges), end(ranges))] <- 1:length(ranges)
	left.is.anchor <- rank[left] > rank[right] 
	ax <- ifelse(left.is.anchor, left, right)
	tx <- ifelse(left.is.anchor, right, left)
	
	# Cleaning up the rest.
	reo <- order(ax, tx)
	final$pairs <- data.frame(anchor.id=ax, target.id=tx)[reo,]
	final$counts <- final$counts[keep,,drop=FALSE][reo,,drop=FALSE]
	final$region <- ranges2
	rownames(final$pairs) <- NULL
	rownames(final$counts) <- NULL
	attributes(final$counts)$dimnames<-NULL
	return(final)
}

###########################################################################################

dir.create("temp-con")
dir1<-"temp-con/1.h5"
dir2<-"temp-con/2.h5"

samecomp <- function(nreads, cuts, ranges, filter=0L, type="any") { 
	simgen(dir1, nreads, chromos)
	simgen(dir2, nreads, chromos)
	
	out <- connectCounts(c(dir1, dir2), fragments=cuts, regions=ranges, filter=filter, type=type)
	ref <- refline(c(dir1, dir2), cuts=cuts, ranges=ranges, filter=filter, type=type)
	if (!identical(ref$pairs$anchor.id, out@anchor.id)) { stop("mismatch in anchor identities") }
	if (!identical(ref$pairs$target.id, out@target.id)) { stop("mismatch in target identities") }
	if (!identical(ref$counts, counts(out))) { stop("mismatch in counts") }
	if (!identical(ref$region, regions(out))) { stop("mismatch in region output") }	

	return(cbind(head(ref$pairs), head(ref$counts)))
}

set.seed(348752)

# Vanilla comparisons involving the same ranges.
current.cuts <- simcuts(chromos)
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=1))
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=2))
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=5))
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=10))

current.cuts <- simcuts(chromos)
samecomp(200, cuts=current.cuts, ranges=simranges(current.cuts, nranges=1))
samecomp(200, cuts=current.cuts, ranges=simranges(current.cuts, nranges=2))
samecomp(200, cuts=current.cuts, ranges=simranges(current.cuts, nranges=5), filter=5L)
samecomp(200, cuts=current.cuts, ranges=simranges(current.cuts, nranges=10), filter=5L)
samecomp(200, cuts=current.cuts, ranges=simranges(current.cuts, nranges=10), type="within")

current.cuts <- simcuts(chromos)
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=1))
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=2), filter=5L)
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=5), filter=5L)
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=10), filter=5L)
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=10), type="within")

current.cuts <- simcuts(chromos)
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=1))
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=2), type="within")
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=5), filter=20L)
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=10), filter=5L)

# Altering the scope of the ranges.
current.cuts <- simcuts(chromos, min=50, max=100, overlap=4)
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300))
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300), type="within")
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300), filter=5)
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=100, min=100, max=300))
samecomp(100, cuts=current.cuts, ranges=simranges(current.cuts, nranges=200, min=100, max=300))

current.cuts <- simcuts(chromos, min=50, max=100, overlap=2)
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300))
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300), type="within")
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300), filter=5)
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=100, min=100, max=300))
samecomp(500, cuts=current.cuts, ranges=simranges(current.cuts, nranges=200, min=100, max=300))
	
current.cuts <- simcuts(chromos, min=50, max=100)
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300))
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300), type="within")
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=50, min=100, max=300), filter=5)
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=100, min=100, max=300))
samecomp(1000, cuts=current.cuts, ranges=simranges(current.cuts, nranges=200, min=100, max=300))

###########################################################################################
# Cleaning up.

unlink("temp-con", recursive=TRUE)

###########################################################################################

