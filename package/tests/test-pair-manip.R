####################################################################################################
# This script tests the pair manipulation functions of savePairs and mergePairs.

suppressPackageStartupMessages(require(diffHic))
suppressPackageStartupMessages(require(rhdf5))

countcomp<-function(n, nfrags, maxc) {
	ai<-as.integer(runif(n, 1, nfrags))
	ti<-as.integer(runif(n, 1, nfrags))
	co<-as.integer(runif(n, 1, maxc))
	out<-data.frame(anchor.id=ai, target.id=ti, count=co)
	out<-out[order(ai, ti),]

	# Checking counting.
	counted<-diffHic:::.sortedAggregate(out)
	comp<-aggregate(co~ai+ti, data=NULL, FUN=length)
	comp<-comp[order(comp[,1], comp[,2]),]
	stopifnot(identical(counted$count, comp[,3]))
	
	# Checking summation.
	counted<-diffHic:::.sortedAggregate(out, mode="sum")
	comp<-aggregate(co~ai+ti, data=NULL, FUN=sum)
	comp<-comp[order(comp[,1], comp[,2]),]
	stopifnot(identical(counted$count, comp[,3]))
	return(head(counted))
}

set.seed(1452312)

countcomp(100, 10, 5)
countcomp(100, 10, 15)
countcomp(100, 10, 25)
countcomp(10, 100, 5)
countcomp(10, 100, 15)
countcomp(10, 100, 25)
countcomp(50, 50, 5)
countcomp(50, 50, 15)
countcomp(50, 50, 25)

####################################################################################################
# Now, checking the behaviour of savePairs. In particular, looking for correct indexing.

tmp<-"temp-pairs"
dir.create(tmp)
savecomp<-function(n, nfrags, nchrs) {
	# Simulating the dudes (target<=anchor at all times).
	ai<-as.integer(runif(n, 1, nfrags))
	ti<-as.integer(runif(n, 1, nfrags))
	out<-data.frame(anchor.id=pmax(ai, ti), target.id=pmin(ai, ti))
	out<-out[order(out$anchor.id, out$target.id),]
	collected<-diffHic:::.sortedAggregate(out)

	# Shuffling the counts.
	original <- collected
	reorder <- nrow(collected):1 
	collected <- collected[reorder,]
	keep <- 1:nrow(collected) %% 3 == 0L
	temp <- collected$anchor.id[keep]
	collected$anchor.id[keep] <- collected$target.id[keep]
	collected$target.id[keep] <- temp

	# Simulating the fragment IDs.
	blah<-GRanges(sample(paste0("chr", 1:nchrs), nfrags, replace=TRUE), IRanges(1:nfrags, 1:nfrags+10),
		seqinfo=Seqinfo(seqnames=paste0("chr", 1:nchrs)))
	blah<-sort(blah)
	newdir<-file.path(tmp, "output")
	savePairs(collected, newdir, fragments=blah)

	# Checking if everything makes sense.
	chrs<-as.character(seqnames(blah))
	indices <- h5ls(newdir)
	indices <- indices[indices$otype=="H5I_DATASET",]
	regot <- list()	
	for (x in 1:nrow(indices)) {
		reread<-h5read(newdir, file.path(indices$group[x], indices$name[x]))
		for (y in 1:ncol(reread)) { attributes(reread[,y]) <- NULL }
		regot[[x]] <- reread
		comp<-diffHic:::.sortedAggregate(reread, mode="sum")
		stopifnot(identical(comp, reread))

		uniq.a<-unique(chrs[reread[,1]])
		uniq.t<-unique(chrs[reread[,2]])
		if (length(uniq.a)!=1L || length(uniq.t)!=1L) { stop("file contains more than one combination") }			
		if (basename(indices$group[x])!=uniq.a || indices$name[x]!=uniq.t) { stop("file contains the incorrect combination") }
	}

	# Checking that the stored result is the same.
	regot <- do.call(rbind, regot)
	regot <- regot[order(regot$anchor.id, regot$target.id),]
	rownames(original) <- rownames(regot) <- NULL 
	stopifnot(identical(original, regot))
	head(regot)
}

savecomp(100, 10, 5)
savecomp(100, 10, 15)
savecomp(100, 10, 25)
savecomp(10, 100, 5)
savecomp(10, 100, 15)
savecomp(10, 100, 25)
savecomp(50, 50, 5)
savecomp(50, 50, 15)
savecomp(50, 50, 25)

####################################################################################################
# Finally, chekcing the merging algorithms.

mergecomp<-function(nl, n, nfrags, nchrs) {
	blah<-GRanges(sample(paste0("chr", 1:nchrs), nfrags, replace=TRUE), IRanges(1:nfrags, 1:nfrags+10),
		seqinfo=Seqinfo(seqnames=paste0("chr", 1:nchrs)))
	blah<-sort(blah)
	allfiles<-list()
	allcounts<-list()
	for (x in 1:nl) {
		# Simulating the dudes (target<=anchor at all times).
		ai<-as.integer(runif(n, 1, nfrags))
		ti<-as.integer(runif(n, 1, nfrags))
		out<-data.frame(anchor.id=pmax(ai, ti), target.id=pmin(ai, ti))
		out<-out[order(out$anchor.id, out$target.id),]
		collected<-diffHic:::.sortedAggregate(out)
		allcounts[[x]]<-collected
		allfiles[[x]]<-  file.path(tmp, paste0("output_", x))
		savePairs(collected, allfiles[[x]], fragments=blah)
	}
	
	# Comparing the combined with a more brutal merger.		
	allfiles<-unlist(allfiles)
	allcounts<-do.call(rbind, allcounts)
	allcounts<-allcounts[order(allcounts[,1], allcounts[,2]),]
	allcounts<-diffHic:::.sortedAggregate(allcounts, mode="sum")
	mdir<-file.path(tmp, "output_merged")
	mergePairs(allfiles, mdir)
	rdir<-file.path(tmp, "output_ref")
	savePairs(allcounts, rdir, fragments=blah)

	# Comparing internal objects.
	combodirs<-c(mdir, rdir)
	out<-diffHic:::.loadIndices(combodirs)
	for (x in names(out)) {
		for (y in names(out[[x]])) {
			current<-out[[x]][[y]]
			stopifnot(all(current))
			stopifnot(identical(h5read(mdir, file.path(x, y)), h5read(rdir, file.path(x, y))))
		}
	}

	return(head(allcounts))
}

mergecomp(2, 100, 10, 5)
mergecomp(3, 100, 10, 15)
mergecomp(4, 100, 10, 25)
mergecomp(3, 10, 100, 5)
mergecomp(4, 10, 100, 15)
mergecomp(2, 10, 100, 25)
mergecomp(4, 50, 50, 5)
mergecomp(2, 50, 50, 15)
mergecomp(3, 50, 50, 25)

####################################################################################################
# Cleaning up.

unlink(tmp, recursive=TRUE)

####################################################################################################


