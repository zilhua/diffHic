scaleAb <- function(data, scale=1, prior.count=2, ...) 
# This computes abundances of `data`, while adjusting for the relative `scale`.
# Mbp. Some care is required regarding the treatment of the prior, here, such
# that the effect of the prior is the same for different scales.
#
# written by Aaron Lun
# 30 October, 2014
{
	if (any(scale <= 0)) { stop("scaling factor should be positive") }
	aveLogCPM(data, prior.count=prior.count*scaled, ...) - log2(scaled)
}
