# Defines the DIList class, that will be used to hold various information.

setClass("DIList", representation(counts="matrix", 
			anchor.id="integer", target.id="integer", region="GRanges",
			coldata="DataFrame", exptdata="List"))

setValidity("DIList", function(object) {
	if (nrow(object@counts)!=length(object@anchor.id)) {
		return('rows in count matrix not equal to length of anchor vector')
	} 
	if (nrow(object@counts)!=length(object@target.id)) { 
		return('rows in count matrix not equal to length of target vector')
	}
	if (ncol(object@counts)!=nrow(object@coldata)) {
		return('columns of count matrix not equal to rows of column data frame')
	}

	if (!all(object@anchor.id >= 1L)) { 
		return('not all anchors are positive integers')
	} 
	if (!all(object@target.id >= 1L)) {
		return('not all targets are positive integers')
	}
	if (!all(object@anchor.id <= length(object@region))) {
		return('not all anchors refer to valid regions')
	} 
	if (!all(object@target.id <= length(object@region))) { 
		return('not all targets refer to valid regions')
	}
	if (!all(object@anchor.id >= object@target.id)) { 
		return('target indices cannot be greater than anchor indices')
	}
	return(TRUE)
})

setMethod("initialize", signature("DIList"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("show", signature("DIList"), function(object) {
	total <- nrow(object@counts)
	nregs <- length(object@region)
	nlibs <- ncol(object@counts)
	cat("DIList object for", nlibs, ifelse(nlibs==1L, "library", "libraries"), 
		"with", total, ifelse(total==1L, "pair", "pairs"), "across", 
		nregs, ifelse(nregs==1L, "region\n", "regions\n"))
})

# Assorted subsetting methods.
setMethod("[", "DIList", function(x, i, j, ..., drop=TRUE) {
	if (missing(i)) {
		new.counts <- x@counts				
		new.anchors <- x@anchor.id
		new.targets <- x@target.id	
	} else {
		new.counts <- x@counts[i,,drop=FALSE]	
		new.anchors <- x@anchor.id[i]
		new.targets <- x@target.id[i]
	}

	if (missing(j)) { 
		new.coldata <- x@coldata
	} else {
		new.counts <- new.counts[,j,drop=FALSE]
		new.coldata <- x@coldata[j,,drop=FALSE]
	}
	initialize(x, counts=new.counts, coldata=new.coldata, 
		anchor.id=new.anchors, target.id=new.targets,
		region=x@region, exptdata=x@exptdata)
})

# Some getters. No need for setters, really.
setGeneric("anchors", function(object, ...) { standardGeneric("anchors") })
setMethod("anchors", signature("DIList"), function(object, id=FALSE) {
	if (id) { return(object@anchor.id) }
	object@region[object@anchor.id]
})

setGeneric("targets", function(object, ...) { standardGeneric("targets") })
setMethod("targets", signature("DIList"), function(object, id=FALSE) {
	if (id) { return(object@target.id) }
	object@region[object@target.id]
})

setGeneric("counts", function(object) { standardGeneric("counts") })
setMethod("counts", signature("DIList"), function(object) {
	object@counts
})

setGeneric("regions", function(object) { standardGeneric("regions") })
setMethod("regions", signature("DIList"), function(object) {
	object@region
})

setMethod("dim", signature("DIList"), function(x) {
	dim(x@counts)
})

setMethod("dimnames", signature("DIList"), function(x) {
	dimnames(x@counts)
})

setMethod("$", signature("DIList"), function(x, name) { 
	x@coldata[[name]]
})

# Borrowing these from GenomicRanges.
setMethod("colData", signature("DIList"), function(x, ...) {
	x@coldata
})

setMethod("exptData", signature("DIList"), function(x, ...) {
	x@exptdata
})

# Constructor object.
DIList <- function(counts, totals=colSums(counts), anchors, targets, regions, expt.data=List(), ...) {
	if (!is.integer(counts)) { storage.mode(counts) <- "integer" }
	anchors <- as.integer(anchors)
	targets <- as.integer(targets)
	totals <- as.integer(totals)
	new("DIList", counts=counts, anchor.id=anchors, target.id=targets, region=regions,
		coldata=DataFrame(totals=totals, ...), exptdata=expt.data)
}

setMethod("c", signature("DIList"), function (x, ..., add.totals=TRUE, recursive=FALSE) {
	if (!identical(recursive, FALSE)) { 
		stop("recursive argument not supported")
	}
	output <- list(counts(x))
	out.a <- list(anchors(x, id=TRUE))
	out.t <- list(targets(x, id=TRUE))
	totality <- x$totals

	ix <- 2L
	for (i in list(...)) {
		if (!is(i, "DIList")) { 
			stop("elements to be concatenated must be DIList objects")
		}
		if (!identical(regions(x), regions(i))) {
			stop("regions should be identical between DIList objects")
		}

		out.a[[ix]] <- anchors(i, id=TRUE)
		out.t[[ix]] <- targets(i, id=TRUE)
		output[[ix]] <- counts(i)

		if (add.totals) { 
			totality <- totality + i$totals 
		} else if (!identical(i$totals, totality)) { 
			warning("totals are not identical between DIList objects")
		}
		ix <- ix + 1L
	}
	
	coldata <- colData(x)
	coldata$totals <- totality
	new("DIList", counts=do.call(rbind, output), 
		anchor.id=unlist(out.a), target.id=unlist(out.t), region=regions(x),
		exptdata=exptData(x), coldata=coldata)
})

# Setting some methods inspired by equivalents in csaw.
setMethod("asDGEList", signature("DIList"), function(object, ...) {
	DGEList(counts(object), lib.size=object$totals, ...)
})

setMethod("normalize", signature("DIList"), function(object, ...) {
	normalizeCounts(counts(object), lib.sizes=object$totals, ...)
})

########################################################################################
# Defining the pairParam class.

setClass("pairParam", representation(fragments="GRanges", restrict="character", discard="GRanges"))

setValidity("pairParam", function(object) {
	# Checking that the fragments are in some order by chromosome name, 
	# and in the correct order by restriction fragment width.		
	if (length(object@fragments)>1L) { 
		if (anyDuplicated(runValue(seqnames(object@fragments)))) { 
			return('restriction fragments should be sorted by chromosome name')	
		}
		
		unsort <-  diff(start(object@fragments)) <= 0L | diff(end(object@fragments)) <= 0L 
		unsort[head(cumsum(runLength(seqnames(object@fragments))),-1L)] <- FALSE 
		# Should be +1, to get to the first element of each chromosome; but, unsort 
		# is missing the first element (because of diff), so no need to add 1.

		if (any(unsort)) {
			return('restriction fragments should be sorted by start and end coordinates')
		}
	}

	if (any(strand(object@fragments)!="*") ) {
		return('restriction fragment ranges should be unstranded')
	}
	return(TRUE)
})

setMethod("initialize", signature("pairParam"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("$", signature("pairParam"), function(x, name) { 
	slot(x, name)
})

setMethod("show", signature("pairParam"), function(object) {
#	if (is.na(object@min.inward)) { 
#		cat("No minimum insert size specified for inward-facing read pairs\n")
#	} else {
#		cat("Minimum insert size for inward-facing read pairs is", object@min.inward, "bp\n")
#	} 
#	if (is.na(object@min.outward)) { 
#		cat("No minimum insert size specified for outward-facing read pairs\n")
#	} else {
#		cat("Minimum insert size for outward-facing read pairs is", object@min.outward, "bp\n")
#	}
#	if (is.na(object@max.frag)) {
#		cat("No maximum fragment size specified\n")
#	} else {
#		cat("Maximum fragment size is", object@max.frag, "bp\n")
#	}

	nfrags <- length(object@fragments)
	nchrs <- length(runValue(seqnames(object@fragments)))
	cat("Genome contains", nfrags, "restriction", ifelse(nfrags==1L, "fragment", "fragments"), 
		"across", nchrs, ifelse(nchrs==1L, "chromosome\n", "chromosomes\n"))

	ndisc <- length(object@discard)
	if (!ndisc) { 
		cat("No discard regions are specified\n")
	} else {
		cat(ndisc, ifelse(ndisc==1L, "region", "regions"), "specified in which alignments are discarded\n")
	}

	nr <- length(object@restrict)
	if (!nr) { 
		cat("No limits on chromosomes for read extraction\n")
	} else {
		if (!attributes(object@restrict)$only.pair) {
			cat("Read extraction is limited to", nr, ifelse(nr==1L, "chromosome\n", "chromosomes\n"))
		} else {
			cat("Read extraction is limited to pairs between", 
				paste0("'", object@restrict[1], "'"), "and", paste0("'", object@restrict[2], "'\n"))
		}
	}
})

pairParam <- function(fragments, 
#	min.inward=NA, min.outward=NA, max.frag=NA, 
	discard=GRanges(), restrict=NULL)
# This creates a SimpleList of parameter objects, specifying
# how reads should be extracted from the BAM files. The aim is
# to synchronize read loading throughout the package, such that
# you don't have to manually respecify them in each function.
#
# written by Aaron Lun
# 1 September 2014
{
#	max.frag <- as.integer(max.frag)
#	min.inward <- as.integer(min.inward)
#	min.outward <- as.integer(min.outward)
	restrict <- .editRestrict(restrict) 
	new("pairParam", 
#			max.frag=max.frag, min.inward=min.inward, min.outward=min.outward,
		restrict=restrict, discard=discard, fragments=fragments)
}

.editRestrict <- function(restrict) {
	only.pair <- FALSE
	if (!is.null(dim(restrict))) { 
		if (nrow(restrict)!=1L || ncol(restrict)!=2L) {
			stop("restrict matrix can only have a single row with two values")
		}
		only.pair <- TRUE
	}
	restrict <- as.character(restrict)
	attr(restrict, "only.pair") <- only.pair
	restrict
}

setGeneric("reform", function(x, ...) { standardGeneric("reform") })
setMethod("reform", signature("pairParam"), function(x, ...) {
	incoming <- list(...)
	sn <- slotNames(x)
	for (sx in names(incoming)) {
		val <- incoming[[sx]]
		sx <- match.arg(sx, sn)
		incoming[[sx]] <- switch(sx, 
#			max.frag=as.integer(val),
#			min.inward=as.integer(val),
#			min.outward=as.integer(val),
			restrict=.editRestrict(val),
			val)
	}
	do.call(initialize, c(x, incoming))
}) 

