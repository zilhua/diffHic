savePairs <- function(x, file, fragments)
# This function saves all counts in 'x' into a set of gzipped files in 'dir', along with
# an index file specifying the identity of each observed chromosome combination corresponding
# to each file. This speeds up any attempt at random access. The idea is to act as a 
# convenience function if you have a matrix of counts (or whatever) and you just want to save it.
#
# written by Aaron Lun  
{
	swap <- x$anchor.id < x$target.id
	if (any(swap)) { 
		temp <- x$target.id[swap]
		x$target.id[swap] <- x$anchor.id[swap]
		x$anchor.id[swap] <- temp
	}
    delta.a <- diff(x$anchor.id)
    delta.t <- diff(x$target.id)
    if (any(delta.a<0L) || any(delta.a==0L & delta.t<0L)) { 
		x <- x[order(x$anchor.id, x$target.id),]
		delta.a <- diff(x$anchor.id)
		delta.t <- diff(x$target.id)
	}
    if (file.exists(file)) { unlink(file, recursive=TRUE) }

    # Need to reorder so fragments are sorted by chromosome COMBINATION. 
    # Sort is stable, no need to supply x$anchor.id/x$target.id in 'order'.
    frag.out <- .checkFragments(fragments)
	all.chrs <- frag.out$chr
	full.chrs <- rep(1:length(all.chrs), frag.out$end-frag.out$start+1L)
    achr <- full.chrs[x$anchor.id]
    tchr <- full.chrs[x$target.id]
    new.o <- order(achr, tchr)
    x <- x[new.o,]

    # Identifying the lengths of the relevant stretches of chromatin, and saving in assorted directories
    # (one per anchor chromosome).
    out <- .sortedAggregate(data.frame(anchor.id=achr[new.o], target.id=tchr[new.o]))
	.initializeH5(file)
	for (ax in unique(out$anchor.id)) { .addGroup(file, all.chrs[ax]) }
	collected <- list()
    sofar <- 0L
    for (y in 1:length(out$count)) {
        current <- out$count[y]
        .writePairs(x[sofar+1:current,], file, all.chrs[out$anchor.id[y]], all.chrs[out$target.id[y]])
        sofar <- sofar+current
    }
    invisible(NULL)
}
