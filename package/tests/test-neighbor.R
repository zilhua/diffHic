###################################################################################################
# This tests the neighbor-counting code.

suppressPackageStartupMessages(require(diffHic))

comp <- function(npairs, chromos, flanking) {
	flanking <- as.integer(flanking)
	nlibs <- 4L
	lambda <- 5
	nbins <- sum(chromos)
    all.pairs <- rbind(t(combn(nbins, 2)), cbind(1:nbins, 1:nbins))
	all.pairs <- data.frame(anchor.id=all.pairs[,2], target.id=all.pairs[,1])       
   	npairs <- min(npairs, nrow(all.pairs))

	# Setting up some data.
	counts <- do.call(cbind, lapply(1:nlibs, FUN=function(x) { as.integer(rpois(npairs, lambda) + 1) }) )
	data <- list(counts=counts, pairs=all.pairs[sample(nrow(all.pairs), npairs),,drop=FALSE], 
		totals=rep(1, nlibs), region=GRanges(rep(names(chromos), chromos), IRanges(1:nbins, 1:nbins)))
	data$region$nfrags <- rep(1:3, length.out=nbins)
	all.chrs <- as.character(seqnames(data$region))
    last.id <- lapply(split(1:nbins, all.chrs), FUN=max)
    first.id <- lapply(split(1:nbins, all.chrs), FUN=min)

	# Going through every pair of chromosomes.
	out.a <- out.t <- matrix(0, nrow=npairs, ncol=nlibs)
	n.a <- n.t <- integer(npairs)
	for (i in 1:npairs) {
		current.a <- data$pairs$anchor.id[i]
		current.t <- data$pairs$target.id[i]
		upper.a <- pmin(current.a + flanking, last.id[[all.chrs[current.a]]]) 
		lower.a <- pmax(current.a - flanking, first.id[[all.chrs[current.a]]])
		upper.t <- pmin(current.t + flanking, last.id[[all.chrs[current.t]]])
		lower.t <- pmax(current.t - flanking, first.id[[all.chrs[current.t]]])
		
		keep.a <- (data$pairs$target.id==current.t & data$pairs$anchor.id >= lower.a & data$pairs$anchor.id <= upper.a) |
                  (data$pairs$anchor.id==current.t & data$pairs$target.id >= lower.a & data$pairs$target.id <= upper.a)
		keep.t <- (data$pairs$anchor.id==current.a & data$pairs$target.id >= lower.t & data$pairs$target.id <= upper.t) | 
		 		  (data$pairs$target.id==current.a & data$pairs$anchor.id >= lower.t & data$pairs$anchor.id <= upper.t)

		out.a[i,] <- colSums(counts[keep.a,,drop=FALSE])
		out.t[i,] <- colSums(counts[keep.t,,drop=FALSE])
		out.a[i,] <- out.a[i,] - data$counts[i,]
		out.t[i,] <- out.t[i,] - data$counts[i,]
		n.a[i] <- upper.a - lower.a
		n.t[i] <- upper.t - lower.t
	}

	# Converting to integer.
	storage.mode(out.a) <- storage.mode(out.t) <- "integer"
	o <- order(data$pairs$anchor.id, data$pairs$target.id)
	equator <- function(x, y) { length(x)==length(y) && all(abs(x-y) < 1e-8*x) }

	# Checking it over.
	a.counts <- countNeighbors(data, flank=flanking, type="anchor")
	t.counts <- countNeighbors(data, flank=flanking, type="target")
	if (!identical(a.counts$counts, out.a)) { stop("mismatch in counts for anchor background") }
	if (!equator(a.counts$n, n.a)) { stop("mismatch in bin pair numbers for anchor background") }
	if (!identical(t.counts$counts, out.t)) { stop("mismatch in counts for target background") }
	if (!equator(t.counts$n, n.t)) { stop("mismatch in bin pair numbers for target background") }

	b.counts <- countNeighbors(data, flank=flanking, type="both")
	ref.combo <- out.a + out.t
	on.diag <- data$pairs$anchor.id==data$pairs$target.id
	ref.combo[on.diag,] <- ref.combo[on.diag,]/2L
	if (!identical(b.counts$counts, ref.combo)) { stop("mismatch in counts for combined anchor/target background") }
	ref.n <- n.t + n.a
	ref.n[on.diag] <- ref.n[on.diag]/2
	if (!equator(b.counts$n, ref.n)) { stop("mismatch in bin pair numbers for combined anchor/target background") }

	return(head(b.counts$counts))
}

###################################################################################################
# Simulating.

set.seed(3427675)
comp(100, c(chrA=10, chrB=30, chrC=20), 5)
comp(100, c(chrA=10, chrC=20), 5)
comp(100, c(chrA=10, chrB=5, chrC=20), 5)
comp(100, c(chrA=20, chrB=5), 5)

comp(100, c(chrA=10, chrB=30, chrC=20), 10)
comp(100, c(chrA=10, chrC=20), 10)
comp(100, c(chrA=10, chrB=5, chrC=20), 10)
comp(100, c(chrA=20, chrB=10), 10)

comp(200, c(chrA=10, chrB=30, chrC=20), 3)
comp(200, c(chrA=10, chrC=20), 3)
comp(200, c(chrA=10, chrB=5, chrC=20), 3)
comp(200, c(chrA=20, chrB=3), 3)

comp(200, c(chrA=10, chrB=30, chrC=20), 1)
comp(200, c(chrA=10, chrC=20), 1)
comp(200, c(chrA=10, chrB=5, chrC=20), 1)
comp(200, c(chrA=20, chrB=5), 1)

###################################################################################################

