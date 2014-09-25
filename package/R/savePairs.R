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
	new.achr <- achr[new.o]
	new.tchr <- tchr[new.o]
    is.diff <- c(TRUE, diff(new.achr)!=0L | diff(new.tchr)!=0L)
	first.in.combo <- which(is.diff)
	last.in.combo <- c(first.in.combo[-1]-1L, length(new.o))

	.initializeH5(file)
	for (ax in unique(new.achr[first.in.combo])) { .addGroup(file, all.chrs[ax]) }
    for (y in 1:length(first.in.combo)) {
        current <- first.in.combo[y]:last.in.combo[y]
		cur.a <- all.chrs[new.achr[current[1]]] 
		cur.t <- all.chrs[new.tchr[current[1]]]
        .writePairs(x[current,], file, cur.a, cur.t)
    }
    invisible(NULL)
}
