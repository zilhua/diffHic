normalizeCNV <- function(data, margins, ref.col=1, span=0.5, prior.count=3,
 	   split=TRUE, abundance=TRUE)
# This performs two-dimensional loess smoothing, using the counts and the 
# marginal counts to compute the abundance and the marginal fold-changes,
# respectively. Both are used as covariates in the model to smooth out any
# systematic differences in interaction intensity. The aim is to get rid
# of any CNV-induced bias, quantified by the differences in the marginals.
{
	if (abundance) { 
		ab <- aveLogCPM(data$counts, lib.size=data$totals)
	}
	adjc <- log(data$counts + 0.5)
	offsets <- matrix(0, nrow=nrow(data$counts), ncol=ncol(data$counts))
	mab <- cpm(margins$counts, lib.size=margins$totals, log=TRUE, prior.count=prior.count)
	is.matched <- findOverlaps(data$region, margins$region, type="equal", select="first")
	ma.adjc <- mab[is.matched[data$pairs$anchor.id],,drop=FALSE] 
	mt.adjc <- mab[is.matched[data$pairs$target.id],,drop=FALSE]

	for (lib in 1:ncol(data$counts)) {
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
	
		# Fitting a loess curve with the specified covariates.	
		i.fc <- adjc[,lib] - adjc[,ref.col]
		cov.fun <- do.call(lp, all.cov)
		fit <- locfit(i.fc ~ cov.fun, alpha=span, deg=1) 
		offsets[,lib] <- fitted(fit)
	}
	offsets <- offsets - rowMeans(offsets)
	return(offsets)
}
