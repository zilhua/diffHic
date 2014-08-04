###################################################################################################
# This provides some functions to simulation fragments and the count directories.

simgen <- function(dir, num, chromos, max.count=10) {
	bonus<-c(0L, cumsum(as.integer(unlist(chromos))))
	overall<-NULL
	for (i in 1:length(chromos)) { 
		max.anchor<-chromos[[i]];
		suppressWarnings(dir.create(file.path(dir, names(chromos)[i])))
		for (j in 1:i) {
			max.target<-chromos[[j]];
			anchors<-as.integer(floor(runif(num, 1, max.anchor)));
			targets<-as.integer(floor(runif(num, 1, max.target)));
			if (i==j){
				anchor.1<-pmax(anchors, targets);
				target.1<-pmin(anchors, targets);
				anchors<-anchor.1;
				targets<-target.1;
			}
			counts<-as.integer(floor(runif(num, 1, max.count)));
			cyrrebt<-data.frame(anchor.id=anchors+bonus[i], target.id=targets+bonus[j], count=counts); 
			overall<-rbind(overall, cyrrebt)
		}
	}
	tmpfrags<-GRanges(rep(names(chromos), chromos), IRanges(1:sum(chromos), 1:sum(chromos)))
	overall<-overall[order(overall[,1], overall[,2]),]
	overall<-diffHic:::.sortedAggregate(overall, mode="sum")
	savePairs(overall, dir, tmpfrags)
}

# Spawning a new cut site set-up.

simcuts<-function(chromos, min=1000, max=10000, overlap=0L) {
	cuts<-list()
	overlap <- as.integer(overlap)
	for (i in 1:length(chromos)) { 
		frags<-as.integer(runif(chromos[[i]], min, max)) 
		frag.ends<-cumsum(frags) - 0:(chromos[[i]]-1L)*overlap
		frag.starts<-c(1L, frag.ends[-chromos[[i]]]+1L-overlap)
		cur_chr<-names(chromos)[i]
		cuts[[cur_chr]]<-GRanges(cur_chr, IRanges(frag.starts, frag.ends))
		seqlengths(cuts[[cur_chr]]) <- max(frag.ends)
	}
	names(cuts)<-NULL
	suppressWarnings(cuts<-do.call(c, cuts))
	return(cuts)
}

###################################################################################################

