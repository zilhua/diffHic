cutGenome <- function(bs, pattern, overhang=4L) 
# This finds the target cut sites in the genome. It currently only searches the
# sense strand, which is fine because if the patterns is an inverse palindrome.
# Otherwise, there may be some problems as you'd need to search the reverse
# strand as well.
#
# written by Aaron Lun
# a long time ago. 
# last modified 22 March 2014
{
	if (nchar(pattern)%%2L!=0) { stop("recognition site must be even in size") }
	ps <- DNAString(pattern)
	if (reverseComplement(ps)!=ps) { stop("recognition site must be an inverse palindrome") }
	overhang <- as.integer(overhang)
	if (overhang > nchar(pattern) || overhang < 0L || overhang%%2L!=0) { stop("overhang must be a non-negative even integer that is not greater than pattern length") }
	remainder <- (nchar(pattern)-overhang)/2L

	original <- list()
	if (is(bs, "BSgenome")) {
		ref.names <- seqnames(bs)
		gen <- genome(bs)
	} else {
		bs <- readDNAStringSet(bs)
		ref.names <- names(bs)
		gen <- NA
	}	

	for (chr in ref.names) {
       	x <- matchPattern(pattern, bs[[chr]])
		match.start <- start(x)
		if (is.unsorted(match.start)) { match.start <- sort(match.start) }
       	starts <- c(1L, match.start+remainder)
		chrlen <- length(bs[[chr]])
       	ends <- c(match.start+remainder-1L+overhang, chrlen)

		# Eliminating fragments formed by consecutive matches at the start or end.
		# These will be ssDNA, so there's nothing to prime the fill-in. 
		if (remainder==0L) {
			keep <- NULL
			if (ends[1]==nchar(pattern)) { keep <- -1L }
			N <- length(ends)
			if (starts[N]==chrlen-nchar(pattern)+1L) { keep <- c(keep, -N) }

			if (length(keep)) { 
				starts <- starts[keep]
				ends <- ends[keep]
			} 
			if (length(starts)==0L) { 
				starts <- 1L
				ends <- chrlen
			}
		}

		original[[chr]] <- GRanges(chr, IRanges(starts, ends))
	}
		
	names(original) <- NULL
	suppressWarnings(original <- do.call(c, original))
    seqlevels(original) <- ref.names
    suppressWarnings(seqinfo(original) <- Seqinfo(ref.names, seqlengths=seqlengths(bs), genome=gen))
	return(original)
}

# We completely ignore tricky interpretations of consecutive sites.
# For starters, the 'remainder' is so low that the strands are unlikely to stay stuck together until the fill-in step.
# It's also unclear whether cleavage is even possible when the recognition site is at the very end of the DNA fragment.
# This is far too much technical detail to be relevant in non-pathological cases.
# Consideration is only given to the start and end, above, to pass pairParam ordering checks and avoid redundancies.


