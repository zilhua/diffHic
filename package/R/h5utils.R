.loadIndices <- function(y)
# A quick and dirty function for index loading with multiple libraries. This
# produces a list which describes the necessary HDF5 object corresponding to
# each chromosome combination for each library. 
{
	overall <- list()
	ni<-length(y)
	for (ix in 1:ni) {
		current <- h5ls(y[ix])
		keep <- current$otype=="H5I_DATASET"
		all.anchors <- basename(current$group[keep])
		assorted.info <- current$name[keep]

		current <- split(assorted.info, all.anchors)
		for (ac in names(current)) {
			if (is.null(overall[[ac]])) { overall[[ac]]<-list() }
			subcurrent <- current[[ac]]
			for (tc in subcurrent) {
				if (is.null(overall[[ac]][[tc]])) { overall[[ac]][[tc]] <- logical(ni) }
				overall[[ac]][[tc]][ix] <- TRUE
			}
		}
	}
	return(overall)
}

.getPairs <- function(y, anchor, target) { h5read(y, file.path(anchor, target)) }

.initializeH5 <- function(y) {
	if (file.exists(y)) { unlink(y, recursive=TRUE) } 
	if (!h5createFile(y)) { stop(sprintf("failed to create '%s'", y)) }
	return(invisible(NULL))
}

.addGroup <- function(y, anchor) {
	if (!h5createGroup(y, anchor)) { stop("failed to add '%s' group to '%s'", anchor, y) }
	return(invisible(NULL))
}

.writePairs <- function(pairs, y, anchor, target) {
	rownames(pairs) <- NULL
	if (h5write(pairs, y, file.path(anchor, target))) { stop("failed to add tag pair data to '%s'", y) }
	return(invisible(NULL))
}
