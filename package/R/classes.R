# Defines the DIList class, that will be used to hold various information.

setClass("DIList", representation(counts="matrix", totals="integer", anchor.id="integer", target.id="integer", region="GRanges"))

setValidity("DIList", function(object) {
	if (nrow(object@counts)!=length(object@anchor.id)) {
		return('rows in count matrix not equal to length of anchor vector')
	} 
	if (nrow(object@counts)!=length(object@target.id)) { 
		return('rows in count matrix not equal to length of target vector')
	}
	if (ncol(object@counts)!=length(object@totals)) { 
		return('columns of count matrix not equal to length of totals vector')
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
	leftover <- 0
	if (total > 10) {
		toshow <- 5
		leftover <- total - toshow
	} else {
		toshow <- total
	}
	nregs <- length(object@region)
	nlibs <- ncol(object@counts)
	cat("DIList object for", nlibs, ifelse(nlibs==1L, "library", "libraries"), 
		"with", total, ifelse(total==1L, "pair", "pairs"), "across", 
		nregs, ifelse(nregs==1L, "region\n", "regions\n"))
	cat("\n")
	
	cat("Counts:\n")
	print(head(object@counts, toshow))
	if (leftover) { cat("... and", leftover, "more rows\n") }
	cat("\n")
	
	cat("Totals:\n")
	print(object@totals)
	cat("\n")

	a.show <- object@region[head(object@anchor.id, toshow)]
	cat("Anchors:\n")
	print(data.frame(Chr=seqnames(a.show), Start=start(a.show), End=end(a.show)))
	if (leftover) { cat("... and", leftover, "more rows\n") }
	cat("\n")

	cat("Targets:\n")
	t.show <- object@region[head(object@target.id, toshow)]
	print(data.frame(Chr=seqnames(t.show), Start=start(t.show), End=end(t.show)))
	if (leftover) { cat("... and", leftover, "more rows\n") }
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
		new.totals <- x@totals
	} else {
		new.counts <- new.counts[,j,drop=FALSE]
		new.totals <- x@totals[j]
	}
	initialize(x, counts=new.counts, totals=new.totals, 
		anchor.id=new.anchors, target.id=new.targets,
		region=x@region)
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

setGeneric("totals", function(object) { standardGeneric("totals") })
setMethod("totals", signature("DIList"), function(object) { 
	object@totals
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

# Constructor object.
.DIList <- function(counts, totals=colSums(counts), anchors, targets, regions) {
	if (!is.integer(counts)) { storage.mode(counts) <- "integer" }
	anchors <- as.integer(anchors)
	targets <- as.integer(targets)
	totals <- as.integer(totals)
	new("DIList", counts=counts, totals=totals, anchor.id=anchors, target.id=targets, region=regions)
}

# Setting some methods inspired by equivalents in csaw.
setMethod("normalize", signature("DIList"), function(object, ...) {
	normalizeCounts(counts(object), lib.sizes=totals(object), ...)
})

setMethod("average", signature("DIList"), function(object, ...) {
	aveLogCPM(counts(object), lib.size=totals(object), ...)
})

setMethod("asDGEList", signature("DIList"), function(object, ...) {
	DGEList(counts(object), lib.size=totals(object), ...)
})

