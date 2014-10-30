aveMb2 <- function(data, areas=NULL, prior.count=2, ...) 
# This computes abundances of `ref`, standardized to an average bin width of 1
# Mbp. Some care is required regarding the treatment of the prior, here, such
# that the effect of the prior is the same for different observed widths.	
#
# written by Aaron Lun
# 30 October, 2014
{ 
	if (is.null(areas)) { areas <- getArea(data, bp=TRUE) }
	scaled <- areas/1e12
	aveLogCPM(asDGEList(data), prior.count=prior.count*scaled, ...) - log2(scaled)
}
