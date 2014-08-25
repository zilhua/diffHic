preparePairs<-function(bam, fragments, file, dedup=TRUE, yield=1e7, ichim=TRUE, minq=NA)
# This function prepares Hi-C data by stripping out all valid pairs from the BAM file and 
# returning a table describing the interacting fragments of that pair. Diagnnostic data is
# also returned describing various bits and pieces of hiC quality.
#
# written by Aaron Lun
# 30 May, 2013
{
	# Preparing cuts; start positions, end positions, index in 'fragments', segmented by chromosome.
	# Anchor/target order is defined by the order of chromosomes in 'fragments'; earlier chromosomes
	# are designated as targets when compared to later chromosomes.
	scuts <- ecuts <- list()
	boost.idx<- list()
	frag.data <- .checkFragments(fragments)
	chrs <- frag.data$chr

	curends<-end(fragments)
	curstarts<-start(fragments)
	for (x in 1:length(chrs)) {
		curdex <- frag.data$start[x]:frag.data$end[x]
		scuts[[x]] <- curstarts[curdex]
		ecuts[[x]] <- curends[curdex]
		if (x==1L) { 
			boost.idx[[chrs[x]]] <- 0L
		} else {
			boost.idx[[chrs[x]]] <- frag.data$end[x-1L]
		}
	}

	# Checking consistency between SAM chromosome lengths and the ones in the cuts.
    chromosomes<-scanBamHeader(bam)[[1]]$targets 
	if (!all(names(chromosomes) %in% chrs)) { stop("missing chromosomes in cut site list"); }
	for (x in 1:length(chrs)) {
		if (chromosomes[[chrs[x]]]!=tail(ecuts[[x]], 1)) {
			stop("length of ", x, " is not consistent between BAM file and fragment ranges")
		}
	}

	# Enforcing input types.
	minq <- as.integer(minq)
	ichim <- as.logical(ichim)
	dedup <- as.logical(dedup)

	# Setting up storage vectors for diagnostics and other output.	
	diagnostics<- 0L
	same.id <- 0L
	singletons <- 0L
	chimeras <- 0L
	allfiles<-list()
	file.count<-1L
	
	# Running through all pairs. Note that, despite the yieldSize, elements are extracted so
	# that runs of 'QNAME's are not broken. See:
	# 	https://stat.ethz.ch/pipermail/bioconductor/2013-March/051490.html for more details.
	bf<-open(BamFile(bam, yieldSize=yield, obeyQname=TRUE, index=character(0)))
	dir <- tempfile(tmpdir=".")
	on.exit(unlink(dir, recursive=TRUE))
	dir.create(dir)

	while (1) {
		out <- scanBam(bf, param=ScanBamParam(what=c("qname", "flag", "rname", "pos", "mapq", "cigar")))[[1]]
		if (!length(out[[1]])) { break; }

		# Converting chromosome ids, using the full set of chromosomes (not just those in this yield cycle).
		rematched <- match(levels(out$rname), chrs)
		if (any(is.na(rematched))) { stop("unrecognised chromosomes in the BAM file") }
		cur.chrs<-rematched[as.integer(out$rname)]-1L

		# Collating read names. We do it here so that any optimization of string comparisons applies 
		# here, rather than having to modify the C++ code from strcmp to pointer comparisons.
		read.pair.len <- rle(out$qname)$length
		collated <- .Call(cxx_report_hic_pairs, scuts, ecuts, read.pair.len, cur.chrs, out$pos,
			out$flag, out$cigar, out$mapq, !ichim, minq, dedup)
		if (is.character(collated)) { stop(collated) }

		# Adding the statistics.
		diagnostics <- diagnostics + collated[[3]]
		same.id <- same.id + collated[[4]]
 	    singletons <- singletons + collated[[5]]	
		chimeras <- chimeras + collated[[6]]
		out <- NULL	

		# Dumping it into a set of temporary files; avoid multiple 'rbind' calls.
		nonempty<-collated[[1]]
		if (!nrow(nonempty)) { next; }
		for (i in 1:nrow(nonempty)) {
			anchor <- chrs[nonempty[i,1]]
			target <- chrs[nonempty[i,2]]
			if (is.null(allfiles[[anchor]])) { allfiles[[anchor]] <- list() }
			
			current.file <- allfiles[[anchor]][[target]]	
			if (is.null(current.file)) { 
				current.file <- file.path(dir, paste0(file.count, ".gz"))
				allfiles[[anchor]][[target]] <- current.file
				file.count <- file.count+1L
			}	
			fout <- gzfile(current.file, open="ab")
			write.table(file=fout, collated[[2]][[i]], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
			close(fout)
		}	
		collated<-NULL
	}
	close(bf)

	# Parsing the directory, pulling out data.frames. We adjust the a/t indices to get the 
 	# full values, and we sort them by anchor/target. We then save them into a HDF5 file.
	.initializeH5(file)
	for (anchor in names(allfiles)) {
		tfiles <- allfiles[[anchor]]
		.addGroup(file, anchor)
		for (target in names(tfiles)) {
			current.file <- tfiles[[target]]
			out <- read.table(current.file, header=FALSE, colClasses="integer")
			colnames(out) <- c("anchor.id", "target.id", "length", "orientation", "gap")
			out$anchor.id <- out$anchor.id+boost.idx[[anchor]]
			out$target.id <- out$target.id+boost.idx[[target]]
			out <- out[order(out$anchor.id, out$target.id),,drop=FALSE]
			rownames(out)<-NULL
			.writePairs(out, file, anchor, target)
		}		
	}	

	# Returning a list of diagnostic components.
    names(diagnostics)<-c("total", "marked", "filtered", "mapped")
 	names(same.id) <- c("dangling", "self.circle")
	names(chimeras) <- c("total", "mapped", "multi", "invalid")
	return(list(pairs=diagnostics, 
				same.id=same.id,
				singles=singletons,
				chimeras=chimeras))
}

getPairData <- function(file, length=TRUE, orientation=TRUE, gap=TRUE)
# This retrieves the fragment sizes, relative orientations and gaps from each directory produced by 
# preparePairs. This is a convenience function which allows people to avoid loading the entire directory 
# in (or manually parsing the said directory).
#
# written by Aaron Lun
{
	# Picking which elements to pull
	if (!(length|orientation|gap)) { 
		stop("must select at least one attribute") 
	}
	topull <- NULL
	if (length) { topull <- c(topull, "length") }
	if (orientation) { topull <- c(topull, "orientation") }
	if (gap)  { topull <- c(topull, "gap") }
	
	# Pulling out data and merging it.
	allfrags <- list()
	ix <- 1L
	allstuff <- .loadIndices(file)
	for (ax in names(allstuff)) {
		current <- allstuff[[ax]] 
		for (tx in names(current)) { 
			allfrags[[ix]] <- .getPairs(file, ax, tx)[,topull,drop=FALSE]
			ix <- ix + 1L
		}
	}
	allfrags <- do.call(rbind, allfrags)

	# Dealing with some loose attributes.
	for (x in 1:ncol(allfrags)) { attributes(allfrags[,x]) <- NULL }
	rownames(allfrags) <- NULL
	return(allfrags)
}

.checkFragments <- function(fragments) 
# Checking incoming fragments for correctness. Checks that the chromosome names
# are in some reasonable order, and checks that the fragment start/ends are all
# sorted (nested fragments are not allowed).
{
	ref.chrs <- as.character(runValue(seqnames(fragments)))
	if (anyDuplicated(ref.chrs)) { stop("fragments should be ordered by chromosome") }
	ref.len <- runLength(seqnames(fragments))
	end.index <- cumsum(ref.len)
	start.index <- end.index - ref.len + 1L

	for (x in 1:length(ref.chrs)) { 
		curf <- fragments[start.index[x]:end.index[x]]
		if (is.unsorted(start(curf)) || is.unsorted(end(curf))) { 
			stop("fragments per chromosome should be sorted by genomic coordinate") 
		}
	}

	return(list(chr=ref.chrs, start=start.index, end=end.index))
}