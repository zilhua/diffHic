###################################################################################################
# This tests the neighbor-counting code.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

.getLimits <- function(x, flank, start, end) {
	lower.x <- x - flank
	upper.x <- x + flank
	if (lower.x < start) { upper.x <- upper.x + start - lower.x }
	if (upper.x > end) { lower.x <- lower.x + end - upper.x }
	lower.x <- max(start, lower.x)
	upper.x <- min(upper.x, end)
	return(c(lower.x, upper.x))
}

comp <- function(npairs, chromos, flanking) {
	flanking <- as.integer(flanking)
	full.flank <- flanking*2L

	nlibs <- 4L
	lambda <- 5
	nbins <- sum(chromos)
    all.pairs <- rbind(t(combn(nbins, 2)), cbind(1:nbins, 1:nbins))
	aid <- pmax(all.pairs[,1], all.pairs[,2])
	tid <- pmin(all.pairs[,1], all.pairs[,2])
   	npairs <- min(npairs, nrow(all.pairs))

	# Setting up some data.
	counts <- do.call(cbind, lapply(1:nlibs, FUN=function(x) { as.integer(rpois(npairs, lambda) + 1) }) )
	chosen <- sample(nrow(all.pairs), npairs)
	indices <- unlist(sapply(chromos, FUN=function(x) { 1:x }), use.names=FALSE)
	data <- DIList(counts=counts, anchors=aid[chosen], targets=tid[chosen],
		totals=rep(1, nlibs), regions=GRanges(rep(names(chromos), chromos), IRanges(indices, indices)))
	data@region$nfrags <- rep(1:3, length.out=nbins)
	
	all.chrs <- as.character(seqnames(regions(data)))
    last.id <- lapply(split(1:nbins, all.chrs), FUN=max)
    first.id <- lapply(split(1:nbins, all.chrs), FUN=min)
	all.diags <- data@anchor.id - data@target.id

	# Going through every pair of chromosomes.
	ref.counts <- matrix(0L, nrow=npairs, ncol=nlibs)
	ref.n <- integer(npairs)
	for (i in 1:npairs) {
		current.a <- data@anchor.id[i]
		current.t <- data@target.id[i]

		if (all.chrs[current.a]==all.chrs[current.t]) {
			chr.start <- first.id[[all.chrs[current.t]]] 
			chr.end <- last.id[[all.chrs[current.t]]] 

			# Diagonal-based background for the intra-chromosomals.
			cur.diag <- all.diags[i]
			t.out <- .getLimits(current.t, flanking, chr.start, chr.end - cur.diag)
			lower.t <- t.out[1]
			upper.t <- t.out[2]
	
			keep <- (all.diags == cur.diag & data@target.id >= lower.t & data@target.id <= upper.t & data@target.id!=current.t)
			ref.counts[i,] <- as.integer(colSums(counts[keep,,drop=FALSE])+0.5)
 			ref.n[i] <- upper.t - lower.t

			# Implementing remedial action. 
			difference <- full.flank - ref.n[i]
			if (difference > 0L) {
				extra.left <- extra.right <- as.integer(difference/2L)
				if (difference %% 2L == 1L) {
					if (current.t <= (chr.start+chr.end-cur.diag+1L)/2L) {
						extra.left <- extra.left + 1L					
					} else {
						extra.right <- extra.right + 1L
					}
				}
		
				# Adding boxes on the edges.
				keep <- (all.diags <= cur.diag - 1L & all.diags >= cur.diag - extra.left & data@target.id == chr.start)
				if (any(keep)) { ref.counts[i,] <- ref.counts[i,] + as.integer(colSums(counts[keep,,drop=FALSE])+0.5) } 
#				print(paste("My own stats are", cur.diag, current.t-chr.start))
#				if (any(keep)){ 
#					print("Left adding:")
#					print(all.diags[keep])
#					print(data@target.id[keep]-chr.start)
#				}
				
				ref.n[i] <- ref.n[i] + ifelse(extra.left > cur.diag, cur.diag, extra.left)
				keep <- (data@anchor.id == chr.end & all.diags <= cur.diag - 1L & all.diags >= cur.diag - extra.right)
				if (any(keep)) { ref.counts[i,] <- ref.counts[i,] + as.integer(colSums(counts[keep,,drop=FALSE])+0.5) }
#				if (any(keep)) {
#					print("Right adding:")
#					print(all.diags[keep])
#					print(data@target.id[keep]-chr.start)
#				}
				ref.n[i] <- ref.n[i] + ifelse(extra.right > cur.diag, cur.diag, extra.right)
			}
		} else {
			# Cross-based background for the inter-chromosomals.
			a.out <- .getLimits(current.a, flanking, first.id[[all.chrs[current.a]]], last.id[[all.chrs[current.a]]])
			lower.a <- a.out[1]
			upper.a <- a.out[2]
			t.out <- .getLimits(current.t, flanking, first.id[[all.chrs[current.t]]], last.id[[all.chrs[current.t]]])
			lower.t <- t.out[1]
			upper.t <- t.out[2]

			keep <- (data@target.id==current.t & data@anchor.id >= lower.a & data@anchor.id <= upper.a & data@anchor.id!=current.a) |
				    (data@anchor.id==current.a & data@target.id >= lower.t & data@target.id <= upper.t & data@target.id!=current.t)  
			ref.counts[i,] <- as.integer(colSums(counts[keep,,drop=FALSE])+0.5) 
			ref.n[i] <- upper.a - lower.a + upper.t - lower.t
		}
	}

	# Converting to integer.
	bg <- countNeighbors(data, flank=flanking)
	obs <- as.matrix(bg$y$counts)
	dimnames(obs) <- NULL
#	print(data[which(obs[,1]!=ref.counts[,1]),])
#	print(obs[which(obs[,1]!=ref.counts[,1]),])
#	print(ref.counts[which(obs[,1]!=ref.counts[,1]),])
	if (!identical(obs, ref.counts)) { stop("mismatch in counts for combined anchor/target background") }
	if (!identical(bg$n, ref.n)) { stop("mismatch in bin pair numbers for combined anchor/target background") }
#	print("YAY")
	return(head(obs))
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

