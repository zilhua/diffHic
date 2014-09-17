###################################################################################################
# This tests the interaction counting capabilities of, well, the interaction counter. 

set.seed(1001)
Sys.setlocale(category="LC_COLLATE",locale="C")
chromos<-c(chrA=51, chrB=31)
source("simcounts.R")

binid <- function(curcuts, dist) {
	# Need to figure out which fragments belong in which bins.
	mids<-(start(curcuts)+end(curcuts))/2
	as.integer((mids-0.1)/dist)+1L
}

refnames<-c("count1", "count2", "anchor", "target")

# We set up the comparison function to check our results. 

finder <- function(dir1, dir2, dist, cuts, filter=10L, restrict=NULL) {
	overall<-list()
	odex<-1L
	totals<-c(0L, 0L)
	collected<-c(0, cumsum(chromos))

	# We need to determine who's who.
	x1 <- h5ls(dir1)
	x1 <- x1[x1$otype=="H5I_DATASET",]
	x1 <- data.frame(anchor=basename(x1$group), target=x1$name)
	x2 <- h5ls(dir2)
	x2 <- x2[x2$otype=="H5I_DATASET",]
	x2 <- data.frame(anchor=basename(x2$group), target=x2$name)
	allbins<-binid(cuts, dist)

	for (k in 1:length(chromos)) {
		cur.k<-names(chromos)[k]
		if (!is.null(restrict) && !cur.k %in% restrict) { next }
		current.krange<-cuts[cur.k==seqnames(cuts)]
		kbin<-binid(current.krange, dist)
		krle <- rle(kbin)
		kend <- cumsum(krle$length)
		kstart <- kend - krle$length + 1L

		for (l in 1:k) {
			cur.l<-names(chromos)[l]
			if (!is.null(restrict) && !cur.l %in% restrict) { next }
			current.lrange<-cuts[seqnames(cuts)==cur.l]
			lbin<-binid(current.lrange, dist)

			# Loading counts.
			x <- list()
			for (m in 1:2) { 
				x.current <- get(paste0("x", m))
				if (any(x1$anchor==cur.k & x1$target==cur.l)) { 
					x[[m]] <- h5read(get(paste0("dir", m)), file.path(cur.k, cur.l))
				} else {
					x[[m]] <- data.frame(anchor.id=integer(0),target.id=integer(0), count=integer(0))
				}
			}

			max.anchor<-chromos[[k]]
			anti.a<-collected[k]
			max.target<-chromos[[l]]
			anti.t<-collected[l]
			mat1<-mat2<-matrix(0L, nrow=max.anchor, ncol=max.target);
			for (g in 1:length(x)) {
				for (i in 1:nrow(x[[g]])) {
					a<-x[[g]][i,1]-anti.a
					t<-x[[g]][i,2]-anti.t
					m<-x[[g]][i,3]
					if (g==1) { 
						mat1[a,t]<-mat1[a,t]+m 
						totals[1]<-totals[1]+m
					} else { 
						mat2[a,t]<-mat2[a,t]+m 
						totals[2]<-totals[2]+m
					}
				}
			}

			# Computing the region for each set of bins (assuming cuts are sorted).
			astart<-start(current.krange)[kstart]
			aend<-end(current.krange)[kend]
			lrle <- rle(lbin)
			lend <- cumsum(lrle$length)
			lstart <- lend - lrle$length + 1L
			tstart <- start(current.lrange)[lstart]
			tend <- end(current.lrange)[lend]
			
			# Computing the counts for each set of bins.
			counts1<-split(mat1, kbin)
			counts2<-split(mat2, kbin)
			current.counts<-list()
			current.anchor<-current.target<-list()
			current.area <- list()
			idex<-1L

			for (g in 1:length(counts1)) {
				subcounts1<-split(matrix(counts1[[g]], nrow=max.target, byrow=TRUE), lbin)
				subcounts2<-split(matrix(counts2[[g]], nrow=max.target, byrow=TRUE), lbin)
				for (i in 1:length(subcounts1)) {
					count.pair<-c(sum(subcounts1[[i]]), sum(subcounts2[[i]]))
					if (sum(count.pair)<filter) { next }
					current.counts[[idex]]<-count.pair
					current.anchor[[idex]]<-c(astart[g], aend[g])
					current.target[[idex]]<-c(tstart[i], tend[i])
					if (k==l && g==i) { 
						current.area[[idex]] <- krle$length[g] * (krle$length[g] + 1)/2
					} else {
						current.area[[idex]] <- krle$length[g] * lrle$length[i]
					}
					idex<-idex+1L
				}
			}
		
			# Aggregating all data.
			if (length(current.counts)) {
				tempa<-do.call(rbind, current.anchor)
				tempt<-do.call(rbind, current.target)
				out<-data.frame(do.call(rbind, current.counts), paste0(cur.k, ":", tempa[,1], "-", tempa[,2]),
						paste0(cur.l, ":", tempt[,1], "-", tempt[,2]), stringsAsFactors=FALSE)
				overall[[odex]]<-out
				odex<-odex+1L
			}
		}
	}

	if (length(overall)) { 
		overall<-do.call(rbind, overall)
		overall<-overall[do.call(order, overall),]
		rownames(overall)<-NULL
	} else {
		overall<-data.frame(integer(0), integer(0), character(0), 
				character(0), numeric(0), stringsAsFactors=FALSE)
	}
	colnames(overall)<-refnames
	return(list(table=overall, total=totals))
}

suppressPackageStartupMessages(require(diffHic))
suppressPackageStartupMessages(require(rhdf5))

dir.create("temp-inter")
dir1<-"temp-inter/1.h5"
dir2<-"temp-inter/2.h5"

comp<-function(npairs1, npairs2, dist, cuts, filter=1L, restrict=NULL) {
	simgen(dir1, npairs1, chromos)
	simgen(dir2, npairs2, chromos)
	y<-squareCounts(c(dir1, dir2), fragments=cuts, width=dist, filter=filter, restrict=restrict)

	ar <- anchors(y)
	tr <- targets(y)
	if (nrow(y)) {
		overall<-data.frame(counts(y), paste0(as.character(seqnames(ar)), ":", start(ar), "-", end(ar)),
			paste0(as.character(seqnames(tr)), ":", start(tr), "-", end(tr)), stringsAsFactors=FALSE)
	} else {
		overall <- data.frame(integer(0), integer(0), character(0), character(0), numeric(0),
			stringsAsFactors=FALSE)
	}
	names(overall)<-refnames
	overall<-overall[do.call(order, overall),]
	rownames(overall)<-NULL

	ref<-finder(dir1, dir2, dist=dist, cuts=cuts, filter=filter, restrict=restrict)
	if (!identical(totals(y), ref$total) || 
		!identical(totals(y), totalCounts(c(dir1, dir2), fragments=cuts, restrict=restrict))) { 
		stop("mismatches in library sizes") 
	}
	if (!identical(overall, ref$table)) { stop("mismatches in counts or region coordinates") }
	if (filter<=1L && !identical(as.integer(colSums(counts(y))+0.5), totals(y))) { 
		stop("sum of counts from binning should equal totals without filtering") }

#	# Checking the fidelity of the filter.
#	unbridled<-squareCounts(c(dir1, dir2), fragments=cuts, width=dist, self.circles=self, filter=0)
#	keep<-rowSums(unbridled$counts)>=filter
#	if (!identical(y$counts, unbridled$counts[keep,])) { stop("mismatches in counts after filtering") }
#	modded<-unbridled$pairs[keep,]
#	rownames(modded)<-NULL
#	if (!identical(y$pairs, modded)) { stop("mismatches in pair IDs after filtering") }
	return(head(overall))
}

###################################################################################################
# Checking a vanilla count.

comp(10, 20, dist=10000, cuts=simcuts(chromos))
comp(10, 20, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(10, 20, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(10, 20, dist=10000, cuts=simcuts(chromos), filter=5)
comp(10, 20, dist=10000, cuts=simcuts(chromos), filter=20)
comp(10, 20, dist=5000, cuts=simcuts(chromos))
comp(10, 20, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(10, 20, dist=5000, cuts=simcuts(chromos), filter=5)
comp(10, 20, dist=5000, cuts=simcuts(chromos), filter=20)
comp(10, 20, dist=1000, cuts=simcuts(chromos))
comp(10, 20, dist=1000, cuts=simcuts(chromos), filter=5)
comp(10, 20, dist=1000, cuts=simcuts(chromos), filter=20)

# Repeating a couple of times.
comp(10, 10, dist=10000, cuts=simcuts(chromos))
comp(10, 10, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(10, 10, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(10, 10, dist=10000, cuts=simcuts(chromos), filter=5)
comp(10, 10, dist=10000, cuts=simcuts(chromos), filter=20)
comp(10, 10, dist=5000, cuts=simcuts(chromos))
comp(10, 10, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(10, 10, dist=5000, cuts=simcuts(chromos), filter=5)
comp(10, 10, dist=5000, cuts=simcuts(chromos), filter=20)
comp(10, 10, dist=1000, cuts=simcuts(chromos))
comp(10, 10, dist=1000, cuts=simcuts(chromos), filter=5)
comp(10, 10, dist=1000, cuts=simcuts(chromos), filter=20)

comp(10, 20, dist=10000, cuts=simcuts(chromos))
comp(10, 20, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(10, 20, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(10, 20, dist=10000, cuts=simcuts(chromos), filter=5)
comp(10, 20, dist=10000, cuts=simcuts(chromos), filter=20)
comp(10, 20, dist=5000, cuts=simcuts(chromos))
comp(10, 20, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(10, 20, dist=5000, cuts=simcuts(chromos), filter=5)
comp(10, 20, dist=5000, cuts=simcuts(chromos), filter=20)
comp(10, 20, dist=1000, cuts=simcuts(chromos))
comp(10, 20, dist=1000, cuts=simcuts(chromos), filter=5)
comp(10, 20, dist=1000, cuts=simcuts(chromos), filter=20)

###################################################################################################
# Another example, a bit more extreme with more overlaps.

comp(50, 20, dist=10000, cuts=simcuts(chromos))
comp(50, 20, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(50, 20, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(50, 20, dist=10000, cuts=simcuts(chromos), filter=5)
comp(50, 20, dist=10000, cuts=simcuts(chromos), filter=20)
comp(50, 20, dist=5000, cuts=simcuts(chromos))
comp(50, 20, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(50, 20, dist=5000, cuts=simcuts(chromos), filter=5)
comp(50, 20, dist=5000, cuts=simcuts(chromos), filter=20)
comp(50, 20, dist=1000, cuts=simcuts(chromos))
comp(50, 20, dist=1000, cuts=simcuts(chromos), filter=5)
comp(50, 20, dist=1000, cuts=simcuts(chromos), filter=20)

comp(30, 30, dist=10000, cuts=simcuts(chromos))
comp(30, 30, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(30, 30, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(30, 30, dist=10000, cuts=simcuts(chromos), filter=5)
comp(30, 30, dist=10000, cuts=simcuts(chromos), filter=20)
comp(30, 30, dist=5000, cuts=simcuts(chromos))
comp(30, 30, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(30, 30, dist=5000, cuts=simcuts(chromos), filter=5)
comp(30, 30, dist=5000, cuts=simcuts(chromos), filter=20)
comp(30, 30, dist=1000, cuts=simcuts(chromos))
comp(30, 30, dist=1000, cuts=simcuts(chromos), filter=5)
comp(30, 30, dist=1000, cuts=simcuts(chromos), filter=20)

###################################################################################################
# A final example which is the pinnacle of extremity.

comp(200, 100, dist=10000, cuts=simcuts(chromos))
comp(200, 100, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(200, 100, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(200, 100, dist=10000, cuts=simcuts(chromos), filter=5)
comp(200, 100, dist=10000, cuts=simcuts(chromos), filter=20)
comp(200, 100, dist=5000, cuts=simcuts(chromos))
comp(200, 100, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(200, 100, dist=5000, cuts=simcuts(chromos), filter=5)
comp(200, 100, dist=5000, cuts=simcuts(chromos), filter=20)
comp(200, 100, dist=1000, cuts=simcuts(chromos))
comp(200, 100, dist=1000, cuts=simcuts(chromos), filter=5)
comp(200, 100, dist=1000, cuts=simcuts(chromos), filter=20)

comp(50, 200, dist=10000, cuts=simcuts(chromos))
comp(50, 200, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(50, 200, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(50, 200, dist=10000, cuts=simcuts(chromos), filter=5)
comp(50, 200, dist=10000, cuts=simcuts(chromos), filter=20)
comp(50, 200, dist=5000, cuts=simcuts(chromos))
comp(50, 200, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(50, 200, dist=5000, cuts=simcuts(chromos), filter=5)
comp(50, 200, dist=5000, cuts=simcuts(chromos), filter=20)
comp(50, 200, dist=1000, cuts=simcuts(chromos))
comp(50, 200, dist=1000, cuts=simcuts(chromos), filter=5)
comp(50, 200, dist=1000, cuts=simcuts(chromos), filter=20)

# Testing some restriction.
comp(50, 200, dist=10000, cuts=simcuts(chromos), restrict="chrB")
comp(50, 200, dist=10000, cuts=simcuts(chromos, overlap=4), restrict="chrA")
comp(50, 200, dist=10000, cuts=simcuts(chromos, overlap=2), restrict="chrA")
comp(50, 200, dist=10000, cuts=simcuts(chromos), filter=5, restrict="chrB")
comp(50, 200, dist=10000, cuts=simcuts(chromos), filter=20, restrict="chrA")
comp(50, 200, dist=5000, cuts=simcuts(chromos), restrict="chrA")
comp(50, 200, dist=5000, cuts=simcuts(chromos, overlap=2), restrict="chrA")
comp(50, 200, dist=5000, cuts=simcuts(chromos), filter=5, restrict="chrB")
comp(50, 200, dist=5000, cuts=simcuts(chromos), filter=20, restrict="chrA")
comp(50, 200, dist=1000, cuts=simcuts(chromos), restrict="chrA")
comp(50, 200, dist=1000, cuts=simcuts(chromos), filter=5, restrict="chrB")
comp(50, 200, dist=1000, cuts=simcuts(chromos), filter=20, restrict="chrA")

##################################################################################################
# Cleaning up.

unlink("temp-inter", recursive=TRUE)

##################################################################################################
# End.
