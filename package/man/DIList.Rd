\name{DIList-class}
\docType{class}

\alias{DIList-class}
\alias{DIList}

\alias{anchors}
\alias{anchors,DIList-method}
\alias{targets}
\alias{targets,DIList-method}
\alias{regions}
\alias{regions,DIList-method}
\alias{counts}
\alias{counts,DIList-method}

\alias{$,DIList-method}
\alias{colData,DIList-method}
\alias{exptData,DIList-method}

\alias{dim,DIList-method}
\alias{dimnames,DIList-method}

\alias{$<-,DIList-method}
\alias{exptData<-,DIList,SimpleList-method}

\alias{[,DIList,ANY,ANY-method}
\alias{show,DIList-method}
\alias{c,DIList-method}
\alias{as.matrix,DIList-method}

\title{DIList class and methods}

\description{
Using the DIList class to store pairwise genomic interactions.
}

\details{
Each DIList stores counts for pairwise genomic interactions.
Slots are defined as:
\describe{
\item{\code{counts}:}{an integer matrix containing the number of read pairs for an interaction in each library}
\item{\code{colData}:}{a \code{DataFrame} object containing library-specific information in each row, e.g., the total number of read pairs}
\item{\code{exptData}:}{a \code{List} object containing data or parameters pertinent to counting}
\item{\code{anchors}:}{an integer vector specifying the index of the first interacting region}
\item{\code{targets}:}{an integer vector specifying the index of the second interacting region}
\item{\code{regions}:}{a \code{GRanges} object containing the coordinates of all interacting regions}
}
Each row of \code{counts} corresponds to an interaction, while each column corresponds to a library.
Each interaction is parameterized by an \code{anchors} and \code{targets} index, pointing to the anchor and target intervals in \code{regions}, respectively.
The anchor interval is defined as that with the higher genomic start coordinate, compared to that of the target interval.
This means that \code{anchors} is never less than \code{targets} to avoid redundant permutations.
}

\section{Constructor}{
\describe{
	\item{}{\code{DIList(counts, totals=colSums(counts), anchors, targets, regions, exptData=List(), ...)}:
	Returns a DIList object where each entry is used to fill the corresponding slots.
	Type coercion is performed as required.
	Arguments in \code{...} are used as columns in \code{colData}.
	}
}
}

\section{Accessors}{
In the code snippets below, \code{x} is a DIList object.

\describe{
	\item{}{\code{anchors(x, id=FALSE)}:
	Get the \code{GRanges} corresponding to the anchor intervals for each interaction.
	If \code{id=TRUE}, indices to the corresponding intervals in \code{regions} are returned instead.
	}
	
	\item{}{\code{targets(x, id=FALSE)}:
	Get the \code{GRanges} corresponding to the target intervals for each interaction.
	If \code{id=TRUE}, indices to the corresponding intervals in \code{regions} are returned instead.
	}

	\item{}{\code{regions(x)}:
	Get the \code{GRanges} for all interacting regions.
	}

	\item{}{\code{counts(x)}:
	Get the matrix of counts for all interactions in all libraries.
	}

	\item{}{\code{colData(x)}:
	Get the \code{DataFrame} of library-specific information.
	}

	\item{}{\code{exptData(x)}:
	Get the \code{List} of experiment-specific information.
	}

	\item{}{\code{x$name}:
	Get the value of \code{colData(x)$name}.
	}

	\item{}{\code{dim(x)}:
	Get the dimensions, i.e., number of interactions (rows) against number of libraries (columns).
	}

	\item{}{\code{dimnames(x)}:
	Get the dimension names. This returns a list of length 2, where each element is \code{NULL} or a character vector.
	}
}
}

\section{Modifiers}{
In the code snippets below, \code{x} is a DIList object.

\describe{ 
	\item{}{\code{x$name <- value}:
	Assign \code{value} to the \code{name} field of the \code{colData} in \code{x}.
	Primarily intended for modification of the library sizes in \code{totals}, but can also be used to store additional library-specific information.
	}
	\item{}{\code{exptData(x) <- value}:
	Assign a \code{SimpleList} object named \code{value} to the \code{exptData} slot of \code{x}.
	Note that this also works with list accessors, e.g., by assigning to \code{exptData(x)$name}.
	}	
}
}

\section{Subsetting and combining}{
In the code snippets below, \code{x} is a DIList object.

\describe{
    \item{}{\code{x[i, j]}:
  	Get count and coordinate data for all interactions \code{i} in libraries \code{j}.
	Either \code{i} or \code{j} can be missing, in which case all interactions or libraries are returned, respectively.
	}

	\item{}{\code{c(x, ..., add.totals=TRUE)}:
	Merge \code{x} with other \code{DIList} objects in \code{...}, by concatenating the count matrices along with the anchor and target indices.
	By default, the totals will be added across all objects to be merged. 
	This can be turned off by setting \code{add.totals=FALSE}, whereby only the totals of \code{x} are used (a warning will be generated if totals are not identical between objects).
	In all cases, objects to be merged should have the same value in the \code{regions} slot.
	}
}
}

\section{Other methods}{
In the code snippets below, \code{x} is a DIList object.
\describe{
	\item{}{\code{show(x)}:
		Shows a summary of the information in \code{x}.
		Specifically, a string is printed that reports the number of libraries, number of interactions and total number of regions in \code{x}.
	}
	\item{}{\code{as.matrix(x, first=NULL, second=first, fill=NULL, ...)}:
		Converts \code{x} into a contact matrix where rows and columns represent regions on the \code{first} and \code{second} chromosomes, respectively.
		Entries of the matrix correspond to an interaction of \code{x}.
		Each entry is extracted from \code{fill}, where \code{fill} is a numeric vector with one value per interaction in \code{x}.
 	   	By default, the average abundance of each interaction in \code{x} is used as the entry.
		For intra-chromosomal spaces, both sides of the diagonal are filled.
		Multiple chromosomes can also be specified in \code{first} or \code{second} (all chromosomes are used by default). 
		Each row and column is named according to the matching index of \code{regions(x)}.
	}
}
}

\author{Aaron Lun}

\examples{
blah <- DIList(counts=matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), ncol=2), 
    totals=c(10L, 10L), anchors=c(1,2,3,4,5,6), targets=c(1,1,2,2,3,3), 
	regions=GRanges("chrA", IRanges(1:6*10, 1:6*10+9)))
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
regions(blah)

blah$totals
colData(blah)
exptData(blah)

blah$totals <- c(10L, 20L)
exptData(blah)$width <- 5
exptData(blah)

c(blah, blah)
c(blah[1:2,], blah[3:4,], add.totals=FALSE)

as.matrix(blah, anchor="chrA")
as.matrix(blah, fill=runif(nrow(blah)))
}

\seealso{
\code{\link{squareCounts}}
}

\keyword{counting}
