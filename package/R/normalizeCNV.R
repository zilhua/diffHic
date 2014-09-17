normalizeCNV <- function(data, margins, ref.col=1, prior.count=3, split=TRUE, abundance=TRUE, 
	degree=1, span=0.3, maxk=500, ...)
# This performs two-dimensional loess smoothing, using the counts and the 
# marginal counts to compute the abundance and the marginal fold-changes,
# respectively. Both are used as covariates in the model to smooth out any
# systematic differences in interaction intensity. The aim is to get rid
# of any CNV-induced bias, quantified by the differences in the marginals.
{
	# Generating covariates.
	adjc <- log(counts(data) + 0.5)
	mab <- cpm(counts(margins), lib.size=totals(margins), log=TRUE, prior.count=prior.count)
	matched <- matchMargins(data, margins)	
	ma.adjc <- mab[matched$amatch,,drop=FALSE] 
	mt.adjc <- mab[matched$tmatch,,drop=FALSE]
	if (abundance) { 
		ab <- aveLogCPM(counts(data), lib.size=totals(data))
	}

	offsets <- matrix(0, nrow=nrow(data), ncol=ncol(data))
	for (lib in 1:ncol(data)) {
		if (lib==ref.col) { next }
		ma.fc <- ma.adjc[,lib] - ma.adjc[,ref.col]
		mt.fc <- mt.adjc[,lib] - mt.adjc[,ref.col]

		if (split) {
			# Anchor/target distinction is arbitrary, so this coerces otherwise-identical 
			# points into the same part of the covariate space.
			ma.fc2 <- pmax(ma.fc, mt.fc) 
			mt.fc2 <- pmin(ma.fc, mt.fc)
			all.cov <- list(ma.fc2, mt.fc2)
		} else { 
			all.cov <- list(ma.fc + mt.fc) 
		}
		if (abundance) { 
			all.cov[[length(all.cov) + 1]] <- ab 
		}
	
		# Fitting a loess surface with the specified covariates.	
		i.fc <- adjc[,lib] - adjc[,ref.col]
		cov.fun <- do.call(lp, c(all.cov, nn=span, deg=degree))
		fit <- locfit(i.fc ~ cov.fun, maxk=maxk, ..., lfproc=locfit.robust) 
		offsets[,lib] <- fitted(fit)
	}
	offsets <- offsets - rowMeans(offsets)
	return(offsets)
}

matchMargins <- function(data, margins) 
# This function just matches the bin pairs in 'data' to the two indices of
# 'margins' that each bin corresponds to.
#
# written by Aaron Lun
# 17 September 2014	
{
	# Checking to ensure that the regions are the same.
	if (!identical(regions(data), regions(margins))) {
		stop("regions must be the same for bin pair and marginal counts") 
	}
	all.indices <- integer(length(data@region))
	all.indices[margins@anchor.id] <- 1:length(margins@anchor.id)
	amatch <- all.indices[data@anchor.id]
	if (any(amatch==0L)) { stop("non-empty anchor in data that is not in margins") }
	tmatch <- all.indices[data@target.id]
	if (any(tmatch==0L)) { stop("non-empty target in data that is not in margins") }
	
	return(data.frame(amatch, tmatch))
}	
