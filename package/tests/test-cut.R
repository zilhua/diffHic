###################################################################################################
# This script is designed to test the restriction site finder. We use the E.coli genome to do so 
# because it's short and we won't be spending ages putting it together.	

suppressPackageStartupMessages(require("diffHic"))

findRestrictionSites <- function(bs, pattern, ref=NULL) {
	require(Biostrings)
	if (nchar(pattern)%%2L!=0L) { stop("Recognition site must be even in size."); }
	ps<-DNAString(pattern);
	if (reverseComplement(ps)!=ps) { stop("Recognition site must be an inverse palindrome."); }
	remainder<-as.integer(nchar(pattern)/2+0.5);
	
	require(BSgenome)
	out<-list();
	for (chr in seqnames(bs)) {
		# Constructing the effective ranges for extraction.
		x<-matchPattern(pattern, bs[[chr]]);
		cuts<-start(x)+remainder-1L;
		out[[chr]]<-c(cuts, length(bs[[chr]]));
	}
	return(out);
}

###################################################################################################

require(BSgenome.Ecoli.NCBI.20080805)

comp<-function(target, overhang=4L) {
	current<-cutGenome(Ecoli, target, overhang=overhang)
	theoretical<-findRestrictionSites(Ecoli, target)

	# Odds and ends.
	overhang<-as.integer(overhang+0.5)
	half.ohang<-as.integer(overhang/2L+0.5)
	half.pattern<-as.integer(nchar(target)/2L)

	# Checking reprocessed cut sites.
	stopifnot(identical(sort(names(theoretical)), sort(seqlevels(current))))
	for (x in names(theoretical)) {
		n<-length(theoretical[[x]])
 	    cuts1<-end(current[x==seqnames(current)])
		cuts1[-n]<-cuts1[-n]-half.ohang
		stopifnot(identical(cuts1, theoretical[[x]]))
	}
	return(head(current))
}

####################################################################################################
# Restriction site must have a 5' overhang.

comp("GGATCC", 4L) # BamHI

comp("GAATTC", 4L) # EcoRI

comp("AACGTT", 2L) # AclI

comp("GCGC", 2L) # HinP1I

comp("GATC", 4L) # MboI

comp("GCGGCCGC", 4L) # NotI

comp("GGCGCGCC", 4L) # AscI

###################################################################################################
# Also some brief checking for fillGenome, though not too much.

dir.create("temp")
comp<-function(target, overhang=4L, spacer=25L) {
	ofile<-"temp/out.fa"
	out<-fillGenome(Ecoli, target, overhang=overhang, outfile=ofile, spacer=spacer)
	current1<-out$original
	current2<-out$modified
	theoretical<-findRestrictionSites(Ecoli, target)

	# Checking sequence file consistency.
	stuff<-readDNAStringSet(ofile)
	chr.len<-sapply(split(end(current2), as.character(seqnames(current2))), FUN=max)			
	for (x in names(theoretical)) {
		stopifnot(length(stuff[[x]])==chr.len[[x]])
		stopifnot(length(stuff[[x]])==tail(theoretical[[x]], 1)+(length(theoretical[[x]])-1)*(spacer+overhang))
	}

	# Checking the first boundary.
	chosen<-which(sapply(theoretical, FUN=length) > 1L)[1]
	if (is.null(chosen)) { return('no sites'); }
	pos<-theoretical[[chosen]][1]
	chosen<-stuff[[names(theoretical)[chosen]]]

	nlen<-nchar(target)
	space.chosen<-subseq(chosen, pos+overhang/2+1, pos+overhang/2+spacer)
	stopifnot(as.character(space.chosen)==paste(rep("N", spacer), collapse=""))

	left.chosen<-subseq(chosen, pos-nlen/2+1, pos+overhang/2)
	right.chosen<-subseq(chosen, pos+overhang/2+spacer+1, pos+overhang+spacer+nlen/2)
	remainder<-(nlen-overhang)/2
	ref<-paste(substr(target, 1, nlen-remainder), substr(target, remainder+1, nlen), sep="")
	stopifnot(as.character(xscat(left.chosen, right.chosen))==ref)

	return(subseq(chosen, pos-nlen/2+1,  pos+overhang+spacer+nlen/2))
}

comp("AAGCTT", 4L) # HindIII

comp("CTAG", 2L) # BfaI

comp("GGCGCGCC", 4L, spacer=10) # AscI

###################################################################################################

unlink("temp", recursive=TRUE);

###################################################################################################
# End.

