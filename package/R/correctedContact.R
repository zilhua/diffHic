correctedContact <- function(data, iterations=50, exclude.local=1, ignore.low=0.02, winsor.high=0.02, average=TRUE, dispersion=0.05)
# This performs the iterative correction method of Mirny et al. (2012) to
# identify the true contact probability of each patch of the interaction
# space. The idea is to use the true contact probability as a filter
# to identify the top set of patches (assuming that most patches represent
# non-specific ligation). Note that for this to work, you need 'filter=1L' 
# during count loading. There is no alternative to loading everything in,
# because you need to compute the truth by dividing and resumming each count.
#
# written by Aaron Lun
# some time ago	
# modified 4 August 2014
{
	# Checking arguments.
	iterations <- as.integer(iterations)
	if (iterations <= 0L) { stop("number of iterations must be a positive integer") }
	ignore.low <- as.double(ignore.low)
	if (ignore.low >= 1) { stop("proportion of low coverage fragments to ignore should be less than 1") }
  	winsor.high <- as.double(winsor.high)
	if (winsor.high >= 1) { stop("proportion of high coverage interactions to winsorize should be less than 1") }
	exclude.local <- as.integer(exclude.local)
    
	# Setting up.
	stopifnot(all(data$pairs$anchor.id >= data$pairs$target.id))
	is.local <- !is.na(getDistance(data))

    # Computing the average counts, if requested. Otherwise, going through each individual library.
	if (average) { 
   		log.lib <- log(data$totals)
		if (length(log.lib)>1L) {
			ave.counts <- exp(edgeR::mglmOneGroup(data$counts, offset=log.lib - mean(log.lib), dispersion=dispersion))
			ave.counts[is.na(ave.counts)] <- 0
		} else {
			ave.counts <- as.double(data$counts[,1])
		}
		out<-.Call(cxx_iterative_correction, ave.counts, data$pairs$anchor.id, data$pairs$target.id, is.local,
			length(data$region), iterations, exclude.local, ignore.low, winsor.high)
 		if (is.character(out)) { stop(out) }
	} else {
		collected.truth <- collected.bias <- collected.max <- list()
		for (lib in 1:length(data$totals)) { 
			nzero <- data$counts[,lib]!=0L
			per.it <-.Call(cxx_iterative_correction, as.double(data$counts[nzero,lib]), data$pairs$anchor.id[nzero], data$pairs$target.id[nzero], 
				is.local[nzero], length(data$region), iterations, exclude.local, ignore.low, winsor.high)
 			if (is.character(per.it)) { stop(per.it) }		

			full.truth <- rep(NA, length(nzero)) 
			full.truth[nzero] <- per.it[[1]]
			collected.truth[[lib]] <- full.truth
			collected.bias[[lib]] <- per.it[[2]]
			collected.max[[lib]] <- per.it[[3]]
		}
		out <- list(do.call(cbind, collected.truth), do.call(cbind, collected.bias), do.call(cbind, collected.max))
	}
	
	names(out) <- c("truth", "bias", "max")
	return(out)
}

# getBias <- function(bias, index, na.action=min)
# # As its name suggests, gets the bias. Some protection is provided against
# # NA values for the bias. These generally correspond to low-abundance regions 
# # for which iterative correction does not provide stable results. We replace 
# # it with the next-best thing, usually the minimum of the non-NA values.
# {
# 	output <- bias[index]
# 	if (any(is.na(output))) { output[is.na(output)] <- na.action(bias[!is.na(bias)]) }
# 	return(output)
# }

