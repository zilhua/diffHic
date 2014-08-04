###################################################################################################
# This tests the bin square summarization method in diffHic.

suppressPackageStartupMessages(require(diffHic))
source("simcounts.R")

chromos<-c(chrA=50, chrB=100, chrC=80)
comp <- function(reference, widths) {
	cutted <- simcuts(chromos, min=20, max=200, overlap=4)
	collected <- list()
	for (w in widths) {
		output <- list()
		bindata <- diffHic:::.getBinID(cutted, w)
		n <- round(runif(1, 20, 100))
		all.a <- as.integer(runif(n, 1, length(bindata$region)))
		all.t <- as.integer(runif(n, 1, all.a))

		output$pairs <- data.frame(anchor.id=all.a, target.id=all.t)
		output$region <- bindata$region
		oname <- paste0("w", w)
		collected[[oname]] <- output
	}	

	output<- do.call(boxPairs, c(collected, fragments=cutted, reference=reference))
	stopifnot(length(output$indices)==length(widths))
	for (x in 1:length(output$indices)) { 
		curdex <- output$indices[[x]]
		curlist <- collected[[x]]

		# Checking that each bin pair is truly nested within its reported parent.
		parent.a <- output$region[output$pairs$anchor.id[curdex]]
		parent.t <- output$region[output$pairs$target.id[curdex]]
		
		current.a <- curlist$region[curlist$pairs$anchor.id]
		current.t <- curlist$region[curlist$pairs$target.id]
		
		if (! all(start(parent.a) <= start(current.a) & end(parent.a) >= end(current.a)) ) { stop("anchor ranges not nested in parent") }
		if (! all(start(parent.t) <= start(current.t) & end(parent.t) >= end(current.t)) ) { stop("target ranges not nested in parent") }
	}

	return(head(output$pairs))
}

set.seed(74653812)
comp(1000, c(100))
comp(1000, c(100, 500))
comp(1000, c(500))
comp(1000, c(500, 1000))
comp(1000, c(100, 500, 1000))

comp(150, c(10))
comp(150, c(10, 50))
comp(150, c(50))
comp(150, c(50, 150))
comp(150, c(10, 50, 150))

comp(500, c(50))
comp(500, c(250))
comp(500, c(50, 250))

####################################################################################################
# End.
