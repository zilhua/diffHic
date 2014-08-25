connectCounts <- function(files, fragments, regions, filter=1L, type="any")
# This counts the number of connections between specified regions in the genome (i.e. between regions
# in 'anchor' and regions in 'target'). This is designed to make it easier to analyze results in terms
# of genes. Note that everything is rounded up to the nearest outside restriction site (or to the
# nearest inside restriction site, depending).
#
# written by Aaron Lun
{
	nlibs <- length(files)
	if (nlibs==0L) { stop("number of libraries must be positive") } 
	if (!is.integer(filter)) { filter<-as.integer(filter) }

	# Figuring out which regions are anchor or targets.
	fdata <- .checkFragments(fragments)
	matched <- match(as.character(seqnames(regions)), fdata$chr)
	if (any(is.na(matched))) { stop("chromosome present in regions and not in fragments") }
	o <- order(matched, start(regions), end(regions))
	regions <- regions[o]

	# Checking out which regions overlap with each fragment.
	if (any(strand(regions)!="*")) { 
		warning("stranded region ranges have no interpretation, coercing unstrandedness") 
		strand(regions) <- "*"
	}
	if (any(strand(fragments)!="*")) { 
		warning("stranded fragment ranges have no interpretation, coercing unstrandedness") 
		strand(fragments) <- "*"
	}
	olaps <- suppressWarnings(findOverlaps(fragments, regions, type=type))
	by.frag <- .retrieveHits(olaps)

	# Setting up output containers.
    full.sizes <- integer(nlibs)
	out.counts <- list()
	out.left <- list()
	out.right <- list()
	idex<-1L

	chrs <- unique(runValue(seqnames(regions)))
    overall<-.loadIndices(files)
	for (anchor in names(overall)) {
		if (!anchor %in% chrs) { next }
		current<-overall[[anchor]]
			
		for (target in names(current)) {
			if (!target %in% chrs) { next }
            pairs <- .baseHiCParser(current[[target]], files, anchor, target)
            full.sizes <- full.sizes+sapply(1:length(pairs), FUN=function(x) { sum(pairs[[x]]$count) })

			# Extracting counts. Running through the fragments and figuring out what matches where.
			out <- .Call(cxx_count_connect, pairs, by.frag$start, by.frag$end, by.frag$hits, filter)
			if (is.character(out)) { stop(out) }
			out.counts[[idex]] <- out[[3]]
			out.left[[idex]] <- out[[1]]
			out.right[[idex]] <- out[[2]]
			idex <-  idex + 1L
		}
	}

	out.counts <- do.call(rbind, out.counts)
	anchor.id <- unlist(out.left)
	target.id <- unlist(out.right)
	o.all <- order(anchor.id, target.id)

	# Generating a new set of regions.
	new.regs <- .redefineRegions(olaps, fragments, regions)
	new.regs$original <- o
	return(.DIList(counts=out.counts[o.all,,drop=FALSE], totals=full.sizes, 
		anchors=anchor.id[o.all], targets=target.id[o.all], regions=new.regs))
}

.retrieveHits <- function(olaps) { 
	start <- end <- integer(queryLength(olaps))
	ok.frags <- queryHits(olaps)
	is.first <- c(TRUE, diff(ok.frags)!=0L)
	start[ok.frags[is.first]] <- which(is.first)
	end[ok.frags[is.first]] <- c(which(is.first)[-1], length(ok.frags)+1L)
	return(list(start=start, end=end, hits=subjectHits(olaps)))
}

.redefineRegions <- function(olaps, fragments, regions) {
	so <- subjectHits(olaps)
	qo <- queryHits(olaps)
	reo <- order(so, qo)
	so <- so[reo]
	qo <- qo[reo]
	
	s.rle <- rle(so)
	r.fin <- cumsum(s.rle$length)
	r.beg <- r.fin - s.rle$length + 1L
	ranges(regions)[s.rle$value] <- IRanges(start(fragments[qo[r.beg]]), end(fragments[qo[r.fin]]))
	# The preceding step is valid because fragments are sorted and non-nested.

	nfrags <- integer(length(regions))
	nfrags[s.rle$value] <- s.rle$length
	regions$nfrags <- nfrags
	return(regions)		
}