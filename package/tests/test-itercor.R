####################################################################################################
# Tests the iterative correction script.

suppressPackageStartupMessages(require(diffHic))
suppressPackageStartupMessages(require(edgeR))
	
comp<- function(npairs, nfrags, nlibs, lambda=5, dispersion=0.05, winsorize=0.02, discard=0.02, locality=1) {
	all.pairs <- rbind(t(combn(nfrags, 2)), cbind(1:nfrags, 1:nfrags))
	all.pairs <- data.frame(anchor.id=all.pairs[,2], target.id=all.pairs[,1])	
	npairs <- min(npairs, nrow(all.pairs))
	counts <- do.call(data.frame, lapply(1:nlibs, FUN=function(x) { rbeta(npairs, lambda, 1) + 0.5 }) )
	data <- list(counts=counts, pairs=all.pairs[sample(nrow(all.pairs), npairs),,drop=FALSE], 
		totals=rep(1, nlibs), region=GRanges(sort(sample(c("chrA", "chrB", "chrC"), nfrags, replace=TRUE)),
			IRanges(1:nfrags, 1:nfrags)))
	
	# Constructing the values.	
	actual.mat<-matrix(0, nfrags, nfrags)
	ave.count <- exp(mglmOneGroup(counts, offset=numeric(nlibs), dispersion=dispersion))
	for (x in 1:nrow(data$pairs)) {
		a<-data$pairs$anchor.id[x]
		t<-data$pairs$target.id[x]
		actual.mat[a,t]<-ave.count[x]
		if (a!=t) { actual.mat[t,a]<-ave.count[x] }
	}
	# Negating local interations.
	if (locality >= 0L){
 		per.chr <- split(1:nfrags, as.integer(seqnames(data$region)))
		for (curloc in 0:locality) {
			failed <- 0L
			for (chr in per.chr) {
				if (length(chr)<=curloc) { 
					failed <- failed + 1L
					next 
				}
				current <- chr[1:(length(chr) - curloc)]
				actual.mat[(current - 1) * nfrags + current + curloc] <- 0
				actual.mat[(current + curloc - 1) * nfrags + current] <- 0
			}
			if (failed==length(per.chr)) { break }
		}		
	}
	
	# Winsorizing.
	temp.mat <- actual.mat[lower.tri(actual.mat, diag=TRUE)]
	is.nonzero <- temp.mat>1e-6
	winsor.val <- max(temp.mat[is.nonzero][sum(is.nonzero) - rank(temp.mat[is.nonzero], ties="first") + 1L > sum(is.nonzero)*winsorize])
	actual.mat[actual.mat > winsor.val] <- winsor.val

	# Checking those that have a low sum of counts.
	frag.sum <- rowSums(actual.mat) 
	not.empty <- frag.sum > 1e-6
	discard.limit <- min(frag.sum[not.empty][rank(frag.sum[not.empty]) > sum(not.empty)*discard])
	to.discard <- frag.sum < discard.limit
	actual.mat[to.discard,] <- 0
	actual.mat[,to.discard] <- 0

	# Iterative correction.
	bias<-rep(1, nfrags)
	iters <- 50
	for (i in 1:iters) {
		additional<-sqrt(rowSums(actual.mat))
		bias <- bias*additional
		additional[additional==0]<-1
		actual.mat<-t(t(actual.mat/additional)/additional)
	}
	
	# Comparing to the reference implementation. We use a fairly gentle threshold for differences,
	# due to the iterative nature of things (and numerical instability and so forth).
	test <- correctedContact(data, dispersion=dispersion, winsor=winsorize, ignore=discard, 
			iterations=iters, exclude.local=locality)
	if (!identical(is.na(test$bias), to.discard)) { stop("invalid biases do not match up") }
	is.okay <- !is.na(test$bias)
	if (any(abs(test$bias[is.okay]-bias[is.okay]) > 1e-6 * bias[is.okay])) { stop("biases do not match up") }
	return(head(bias))
}

####################################################################################################

set.seed(0)

# Varying the number of fragments.

comp(100, 20, 2, discard=0.1)
comp(100, 30, 2, discard=0.1)
comp(100, 40, 2, discard=0.1)

comp(100, 20, 2, winsor=0.05)
comp(100, 30, 2, winsor=0.1)
comp(100, 40, 2, winsor=0.01)

comp(100, 20, 2, locality=0)
comp(100, 30, 2, locality=2)
comp(100, 40, 2, locality=1)

# Trying with fewer reads.

#debug(comp)
comp(10, 20, 2, discard=0.1)
comp(10, 30, 2, discard=0.1)
comp(10, 40, 2, discard=0.1)

comp(20, 20, 2, winsor=0.05)
comp(20, 30, 2, winsor=0.1)
comp(20, 40, 2, winsor=0.01)

comp(10, 20, 2, locality=0)
comp(10, 30, 2, locality=2)
comp(10, 40, 2, locality=1)

# Trying with fewer libraries.

comp(50, 20, 1, discard=0.1)
comp(50, 30, 1, discard=0.1)
comp(50, 40, 1, discard=0.1)

comp(50, 20, 1, winsor=0.05)
comp(50, 30, 1, winsor=0.1)
comp(50, 40, 1, winsor=0.01)

comp(50, 20, 1, locality=0)
comp(50, 30, 1, locality=2)
comp(50, 40, 1, locality=1)

# Trying with no special attention.
comp(50, 20, 1, discard=0, winsor=0, locality=-1)
comp(50, 50, 1, discard=0, winsor=0, locality=-1)
comp(50, 20, 2, discard=0, winsor=0, locality=-1)
comp(50, 20, 2, discard=0, winsor=0, locality=1000)
comp(100, 20, 2, discard=0, winsor=0, locality=1000)

####################################################################################################
# End.
