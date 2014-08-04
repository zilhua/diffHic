fillGenome <- function(bs, pattern, outfile, overhang=4L, spacer=50L) 
# This function takes a BSgenome object, finds all the restriction sites, and fills them in with 
# a spacer. The idea is to encourage fusion mappers to consider breakpoints at restriction sites
# by also supplying the breakpoints. It also promotes ligation of adjacent fragments as 
# fusions rather than as single alignments with a little inserted region. It then returns a 
# a GenomicRangesList, for the coordinates of the fragments on the original and modified genomes.
#
# written by Aaron Lun
{
	spacer<-as.integer(spacer)
	spaced<-paste(rep("N", spacer), collapse="")
	if (file.exists(outfile)) { unlink(outfile) } # Need to erase it if it already exists, as we're appending.
	
	original <- cutGenome(bs, pattern, overhang=overhang)
	modified<-list()
	maxed<-list()
	for (chr in seqnames(bs)) {
		chosen <- seqnames(original)==chr
		starts <- start(original[chosen])
		ends <- end(original[chosen])

        # Extracting and saving the sequences (as fragments).
        y <- getSeq(bs, names=chr, start=starts, end=ends)
		if (!is(y, "DNAStringSet")) { y <- Biostrings::DNAStringSet(y) }
		filled <- DNAStringSet(paste(y, collapse=spaced))
		names(filled) <- chr
        writeXStringSet(filepath=outfile, filled, append=TRUE)
		
		# Getting the fragmentation coordinates on the filled genome.
		n <- length(y)
		eos <- cumsum(width(y))+1:n*spacer
		ends <- eos-spacer
		starts <- c(1L, eos[-n]+1L)
		modified[[chr]] <- GRanges(chr, IRanges(starts, ends))
		maxed[[chr]] <- ends[n]
    }

	# Returning the two sets of GRanges.
	names(modified)<-NULL
	suppressWarnings(modified<-do.call(c, modified))
    seqlevels(modified)<-seqnames(bs)
	maxed<-unlist(maxed)
    suppressWarnings(seqinfo(modified)<-Seqinfo(names(maxed), maxed))

	return(list(original=original, modified=modified))
}


