#normalizeHiC <- function(data, iter=50, ref=1, trim=0.3, prior.count=0.1)
## This computes offsets for all entries in a Hi-C contact matrix for each library. It requires 
## the same parameters that were fed into squareCounts, so that it will return an offset matrix 
## with the correct dimensions. This applies the iterative correction procedure to each pair of
## libraries until the corrected library sums are equal for each fragment.
#{
#	iter<-as.integer(iter+0.5)
#	ref<-as.integer(ref+0.5)
#	lib.sizes<-data$totals
#    prior.count.scaled <- lib.sizes/mean(lib.sizes) * prior.count
#    lib.sizes <- lib.sizes + 2 * prior.count.scaled
#
#	out<-.Call("R_median_correct", data$counts, data$pairs$anchor, data$pairs$target, 
#		iter, ref, trim, prior.count.scaled, PACKAGE="diffHic")
#	if (is.character(out)) { stop(out) }
#	names(out)<-c("offset", "step")
#	return(out);
#}

matchPairs <- function(data1, data2)
# This computes offsets given a 'data' list and a 'binned' list (constructed using a multiple
# of the width used for 'data'). The idea is to get the 'bins' squares containing the 'data'
# squares of interest, and to use the larger counts of the former to then perform scaling-based
# normalization. This assumes that most squares in 'data' are not DE. 
#
# written by Aaron Lun
# 28 Oct 2013
{
 	matching<-findOverlaps(data1$region, data2$region, type="within", select="first")
	if (any(is.na(matching))) { stop("regions in data1 should be a subinterval of a region in data2") }

	o2<-order(data2$pairs$anchor, data2$pairs$target)
	refined.anchor<-matching[data1$pairs$anchor]
	refined.target<-matching[data1$pairs$target]
	ro<-order(refined.anchor, refined.target)

	out<-.Call("R_match_pairs", ro, refined.anchor, refined.target, 
		o2, data2$pairs$anchor, data2$pairs$target, PACKAGE="diffHic")
	if (is.character(out)) { stop(out) }
	return(out)
}

# Trend correction is difficult because of a number of factors:
#    a) Small counts. This stymies the use of robust empirical methods as it's 
#       not always clear how to fit to discrete points in, for example, the MA 
#       space. Analyses become very sensitive to the choice of prior count when
#       lots of small counts are around (e.g. increased spread of M-values,
#       greater/smaller offsets in quantile). Any non-parametric method will 
#       struggle if it does not account for the discreteness of the count data.
#    b) Dominance of small counts. Related to the first one, but compounds it 
#       because the majority of counts are going to be small. This messes up
#       any span based technique as everything has to account for the small counts
#       floating around. For example, 30 million counts with row sums below 1,
#       and under 300000 with row sums above 20? This basically means that your
#       span must always include low counts, which probably isn't that healthy. 
#    c) Spurious trends upon filtering. One might think that the best way to get
#       rid of small counts is to remove them by row sum filtering. However, removal 
#       of such counts can generate trends in the MA plot or bulk differences in
#       the CDFs, especially for differing library sizes. You can only avoid this
#       by removing low counts after normalization, which does relatively well
#       but still exposes low counts to the normalization algorithm.
#
# Solutions at low counts for a non-parametric method like cyclic loess:
#    a) Filter them out. You have to do this anyway, so it's the most convenient.
#       However, you want to normalize first otherwise you'll end up with 
#       a bunch of spurious trends when offsets are truly non-zero.
#    b) Give them a large enough prior count to suppress the generation of 
#       spurious trends. This reduces the pressure of the filter choice,
#       but can lead to problems when the trend is truly large.
#
# This strategy seems to work well save for the situation where batches occur in 
# both libraries of a group. There, using a low prior helps but that's probably
# a case of compensating errors. Slight improvements in loess over TMM are observed
# though I wouldn't stress it out too much; the key is that it is operable and
# provides satisfactory results.

normalizeHiC <- function (counts, lib.sizes=colSums(counts), dispersion=0.05, ...) 
# This computes offsets given a 'data' object and a prior count. Specifically,
# we compute the normalization factors by fitting a loess curve to the M-values
# (between each library and an 'average' library) using the abundances as the 
# covariate. This is a bit more appropriate than using the A-value i.e. geometric 
# mean. The log-transformation is done with a continuity correction.
#
# written by Aaron Lun
# 12 November 2013
{
	ab <- aveLogCPM(counts, lib.size=lib.sizes, dispersion=dispersion)
	offs <- matrix(0, nrow(counts), ncol(counts), byrow=TRUE)
	for (x in 1:ncol(counts)) { offs[,x]<-loessFit(log(counts[,x]+0.5), ab, ...)$fitted }
	offs<-offs-rowMeans(offs)
	return(offs)
}

# .normalizeHiC <- function (counts, lib.sizes=colSums(counts), dispersion=0.2, iter=3, trim=0.1, ...) 
# # This computes offsets given a 'data' object and a prior count. Specifically,
# # we compute the normalization factors by fitting a loess curve to the observed
# # counts, using the abundances as the covariate. We robustify it by applying
# # biweights on the deviances.
# #
# # written by Aaron Lun
# # 25 November 2013
# {
#     require(locfit)
# 	require(edgeR)
# 	ntag <- nrow(counts)
# 	nlib <- ncol(counts)
# 	offs <- matrix(0, ntag, nlib)
# 	devfun<-deviances.function(dispersion)
# 	ab<-mglmOneGroup(counts, dispersion=dispersion, offset=log(lib.sizes))
# 			
# 	for (i in 1:nlib) {
# 		curweights<-rep(1, ntag)
# 		for (j in 1:iter) { 
# 			out <- locfit.raw(x=ab, y=counts[,i], base=ab+log(lib.sizes[i]), 
# 				weights=curweights, family="qpoisson", link="log", ...)
# 			res <- log(data$counts[,i]+0.5)-log(fitted(out))
#             cmad <- median(abs(res))*6
# 			curweights <- pmax(1 - (res/cmad)^2, 1e-6)^2
# 		}
# 		offs[,i]<-log(fitted(out))
# 	}
# 
# 	offs<-offs-rowMeans(offs)
# 	return(offs)
# }

# This doesn't seem to work as well, as it seems to be quite sensitive to the dispersion 
# used for robustness estimates. Probably with good reason, given that you can get
# serious changes in the deviance with differing dispersions, whereas the effect of fiddling
# with the prior count is generally more predictable.
