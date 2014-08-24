# Defines the DIList class, that will be used to hold various information.

setClass("DIList", representation(counts="matrix", totals="integer", anchor.id="integer", target.id="integer", region="GRanges"))

setValidity("DIList", function(object) {
	if (nrow(object@counts)!=length(object@anchor.id)) {
		return('rows in count matrix not equal to length of anchor vector')
	} 
	if (ncol(object@counts)!=length(object@totals)) { 
		return('columns of count matrix not equal to length of totals vector')
	}
	if (nrow(object@counts)!=length(object@target.id)) { 
		return('rows in count matrix not equal to length of target vector')
	}
	if (!all(object@anchor.id >= 1L)) { 
		return('not all anchors are positive integers')
	} 
	if (!all(object@target.id >= 1L)) {
		return('not all targets are positive integers')
	}
	if (max(object@anchor.id) > length(object@region)) {
		return('not all anchors refer to valid regions')
	} 
	if (max(object@target.id) > length(object@region)) { 
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
setMethod("[", c("DIList", "missing", "numeric", "ANY"), function(x, i, j, ..., drop=TRUE) {
	initialize(x, counts=x@counts[,j,drop=FALSE], totals=x@totals[j],
		anchor.id=x@anchor.id, target.id=x@target.id, 
		region=x@region)
})

setMethod("[", c("DIList", "numeric", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) {
	initialize(x, counts=x@counts[i,,drop=FALSE], totals=x@totals,
		anchor.id=x@anchor.id[i], target.id=x@target.id[i], 
		region=x@region)
})

setMethod("[", c("DIList", "numeric", "numeric", "ANY"), function(x, i, j, ..., drop=TRUE) {
	initialize(x, counts=x@counts[i,j,drop=FALSE], totals=x@totals[j],
		anchor.id=x@anchor.id[i], target.id=x@target.id[i], 
		region=x@region)
})

# Some getters. No need for setters, really.
setGeneric("anchors", function(object){ standardGeneric("anchors") })
setMethod("anchors", signature("DIList"), function(object) {
	object@region[object@anchor.id]
})

setGeneric("targets", function(object){ standardGeneric("targets") })
setMethod("targets", signature("DIList"), function(object) {
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

setGeneric("npairs", function(object) { standardGeneric("npairs") })
setMethod("npairs", signature("DIList"), function(object) {
	nrow(object@counts)
})

setGeneric("nlibs", function(object) { standardGeneric("nlibs") })
setMethod("nlibs", signature("DIList"), function(object) {
	ncol(object@counts)
})	

# Constructor object.
.DIList <- function(counts, totals=colSums(counts), anchors, targets, regions) {
	if (!is.integer(counts)) { storage.mode(counts) <- "integer" }
	anchors <- as.integer(anchors)
	targets <- as.integer(targets)
	totals <- as.integer(totals)
	new("DIList", counts=counts, totals=totals, anchor.id=anchors, target.id=targets, region=regions)
}

# Testing:
# blah <- diffHic:::.DIList(matrix(c(1,1,2,2,3,3,4,4), ncol=2), anchors=c(1,2,3,4), targets=c(1,1,2,2), regions=GRanges("chrA", IRanges(10+1:4, 2+20:23)))
# blah[1,]


