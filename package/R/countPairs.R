.sortedAggregate <- function(x, mode=c("count", "sum"))
# This is an efficient version of 'aggregate(count/length~anchor.id+target.id, data=x, FUN=sum/length)'
# which takes advantage of the fact that x is sorted by anchor/target.id values.
{
	if (!nrow(x)) { 
		return(data.frame(anchor.id=integer(0), target.id=integer(0), count=integer(0)))
	}
	mode <- match.arg(mode)
	delta.a <- diff(x$anchor.id)
	delta.t <- diff(x$target.id)
	if (any(delta.a<0L) | any(delta.a==0L & delta.t<0L)) { stop("anchor/target indices should be sorted") }
	last.one <- which(c(delta.a!=0L | delta.t!=0L, TRUE))
		
	if (mode=="count")  { collected <- last.one }
	else if (mode=="sum") { collected <- cumsum(x$count)[last.one] }
	counts <- collected-c(0L, collected[-length(collected)])
	data.frame(anchor.id=x$anchor.id[last.one], target.id=x$target.id[last.one], count=counts)
}

countPairs <- function(file.in, file.out=file.in, max.length=NA, min.ingap=NA, min.outgap=NA)
# This function collates all valid HiC pairs into a set of counts for each interaction. Basically, it 
# takes the directory output by 'prepareHiCPairs' and converts it into counts for each interaction. 
# If the fragment size is too big then we filter it out. We also filter out pairs which are too
# close together (min.*gap) with the difference between the two specifying inward/outward facing pairs.
#
# writtn by Aaron Lun
{
    # Use a temporary file as a placeholder, in case file.out==file.in.
	tmpf <- tempfile(tmpdir=".")
	on.exit({ if (file.exists(tmpf)) { unlink(tmpf) } })
	.initializeH5(tmpf)
	retained <- total <- by.len <- by.in <- by.out <- 0L

	# Parsing through the old index, counting/summing everything, and saving it to the
	# temporary file. We also remove any specified elements.
	allstuff <- .loadIndices(file.in)
	for (ax in names(allstuff)) { 
		current <- allstuff[[ax]]
		loaded <- FALSE
		for (tx in names(current)) { 
			x <- .getPairs(file.in, ax, tx)

			if (!is.na(max.length)) { 
				keep.len <- x$length <= max.length
				by.len <- by.len + sum(!keep.len)
			} else { keep.len <- TRUE }
			if (!is.na(min.ingap)) { 
				keep.in <- x$orientation!=1L | is.na(x$gap) | x$gap >= min.ingap
				by.in <- by.in + sum(!keep.in)
			} else { keep.in <- TRUE }
			if (!is.na(min.outgap)) { 
				keep.out <- x$orientation!=2L | is.na(x$gap) | x$gap >= min.outgap
				by.out <- by.out + sum(!keep.out)
			} else { keep.out <- TRUE }

			total <- total + nrow(x) 
			x <- x[keep.len & keep.in & keep.out,]
			if (nrow(x)) { 
				if (!loaded) { # Only adding a group if the data.frame is non-empty.
					.addGroup(tmpf, ax)
					loaded <- TRUE
				}
				.writePairs(.sortedAggregate(x, mode="count"), tmpf, ax, tx)
				retained <- retained + nrow(x)
			}
		}
	}
	
	# Shuffling things around.
	if (!file.rename(tmpf, file.out)) { stop("cannot move file to the specified destination") }
	return(c(total=total, length=by.len, ingap=by.in, outgap=by.out, retained=retained))
}
