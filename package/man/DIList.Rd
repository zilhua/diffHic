\name{DIList-class}
\docType{class}

\alias{DIList-class}

\alias{anchors}
\alias{anchors,DIList-method}
\alias{targets}
\alias{targets,DIList-method}
\alias{regions}
\alias{regions,DIList-method}
\alias{counts}
\alias{counts,DIList-method}
\alias{totals}
\alias{totals,DIList-method}

\alias{dim,DIList-method}
\alias{dimnames,DIList-method}

\alias{[,DIList,ANY,ANY-method}
\alias{show,DIList-method}

\alias{average}
\alias{average,DIList-method}
\alias{normalize}
\alias{normalize,DIList-method}
\alias{asDGEList}
\alias{asDGEList,DIList-method}

\title{DIList class and methods}

\description{
Using the DIList class to store pairwise genomic interactions.
}

\details{
Each DIList stores counts for pairwise genomic interactions.
Slots are defined as:
\describe{
\item{\code{counts}}{an integer matrix containing the number of read pairs for an interaction in each library}
\item{\code{totals}}{an integer vector containing the total number of read pairs for each library}
\item{\code{anchor.id}}{an integer vector specifying the index of the first interacting region}
\item{\code{target.id}}{an integer vector specifying the index of the second interacting region}
\item{\code{region}}{a \code{GRanges} object containing the coordinates of all interacting regions}
}
Each row of \code{counts} corresponds to an interaction, while each column corresponds to a library.
Each interaction is parameterized by an \code{anchor.id} and \code{target.id} index, pointing to the anchor and target intervals in \code{region}, respectively.
The anchor interval is defined as that with the higher genomic coordinate compared to the target interval, to avoid redundancy from permutations.
}

\section{Accessors}{
In the code snippets below, \code{x} is a DIList object.

\describe{
	\item{}{\code{anchors(x, id=FALSE)}:
	Get the \code{GRanges} corresponding to the anchor intervals for each interaction.
	If \code{id=TRUE}, indices to the corresponding intervals in \code{region} are returned instead.
	}
	
	\item{}{\code{targets(x, id=FALSE)}:
	Get the \code{GRanges} corresponding to the target intervals for each interaction.
	If \code{id=TRUE}, indices to the corresponding intervals in \code{region} are returned instead.
	}

	\item{}{\code{regions(x)}:
	Get the \code{GRanges} for all interacting regions.
	}

	\item{}{\code{counts(x)}:
	Get the matrix of counts for all interactions in all libraries.
	}

	\item{}{\code{totals(x)}:
	Get the vector of total library sizes for all libraries.
	}

	\item{}{\code{dim(x)}:
	Get the dimensions, i.e., number of interactions (rows) against number of libraries (columns).
	}

	\item{}{\code{dimnames(x)}:
	Get the dimension names. This returns a list of length 2, where each element is \code{NULL} or a character vector.
	}
}
}

\section{Subsetting}{
In the code snippets below, \code{x} is a DIList object.

\describe{
    \item{}{\code{x[i, j]}:
  	Get data for all interactions \code{i} in libraries \code{j}.
	Either \code{i} or \code{j} can be missing, in which case all interactions or libraries are returned, respectively.
	}
}
}

\section{Other methods}{
\describe{
	\item{}{\code{show(x)}:
		By default the \code{show} method displays count and interval data for the top 5 interactions.
		Data for all interactions will be shown if there are no more than 10 interactions.
	}
	\item{}{\code{normalize(x, ...)}:
		Normalizes the counts between libraries in \code{x}.
 		Count data and library sizes are passed to the \code{\link{normalizeCounts}} function, along with any arguments in \code{...}.
		A numeric matrix or vector is then returned, depending on the type of normalization.
	}
	\item{}{\code{average(x, ...)}:
		Computes the average log-CPM of each interaction in \code{x}.
		Count data, library sizes and arguments in \code{...} are passed to the \code{\link{aveLogCPM}} function, and a numeric vector is returned.
	}
	\item{}{\code{asDGEList(x, ...)}:
		Extracts count data and library sizes and returns the same data as a \code{DGEList}.
		Arguments in \code{...} are passed to the \code{\link{DGEList}} constructor.
	}
}
}

\author{Aaron Lun}

\examples{
blah <- new("DIList", counts=matrix(as.integer(c(1,1,2,2,3,3,4,4)), ncol=2), totals=c(10L, 10L), 
	anchor.id=as.integer(c(1,2,3,4)), target.id=as.integer(c(1,1,2,2)), 
	region=GRanges("chrA", IRanges(10+1:20, 2+21:40)))
nrow(blah)
ncol(blah)

blah
blah[1,]
blah[,2]
blah[1:2,2]

anchors(blah)
anchors(blah, id=TRUE)
targets(blah)
targets(blah, id=TRUE)
counts(blah)
totals(blah)
regions(blah)

normalize(blah, type="loess")
average(blah)
asDGEList(blah)
}

\seealso{
\code{\link{normalizeCounts}},
\code{\link{aveLogCPM}},
\code{\link{DGEList}}
}
