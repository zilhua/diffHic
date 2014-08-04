\name{cutGenome}
\alias{cutGenome}

\title{Cut up the genome}

\description{Perform an in silico restriction digest of a target genome.}

\usage{
cutGenome(bs, pattern, overhang=4L)
}

\arguments{
	\item{bs}{a \code{BSgenome} object}
	\item{pattern}{character string describing the recognition site}
	\item{overhang}{integer scalar specifying the length of the 5' overhang}
}

\value{
A \code{GRanges} object containing the boundaries of each restriction fragment in the genome.
}

\details{
This function simulates a restriction digestion given the recognition site and
overhang. The total sequence spanned by each fragment is recorded, including
the two sticky ends. No support is currently provided for searching the reverse
strand, so the recognition site should be an inverse palindrome.
}

\author{Aaron Lun}

\examples{
require(BSgenome.Ecoli.NCBI.20080805)

# HindIII
cutGenome(Ecoli, "AAGCTT", overhang=4L)

# SacII
cutGenome(Ecoli, "CCGCGG", overhang=2L)

# AluI
cutGenome(Ecoli, "AGCT", overhang=0L)
}

\seealso{
\code{\link{fillGenome}}
}

