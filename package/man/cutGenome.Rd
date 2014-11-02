\name{cutGenome}
\alias{cutGenome}

\title{Cut up the genome}

\description{Perform an in silico restriction digest of a target genome.}

\usage{
cutGenome(bs, pattern, overhang=4L)
}

\arguments{
	\item{bs}{a \code{BSgenome} object or a character string pointing to a FASTA file}
	\item{pattern}{character string describing the recognition site}
	\item{overhang}{integer scalar specifying the length of the 5' overhang}
}

\value{
A \code{GRanges} object containing the boundaries of each restriction fragment in the genome.
}

\details{
This function simulates a restriction digestion of a specified genome, given the recognition site and 5' overhang of the cutter
The total sequence spanned by each fragment is recorded, including the two sticky ends. 
No support is currently provided for searching the reverse strand, so the recognition site should be an inverse palindrome.

The genome should be specified as a \code{BSgenome} object.
However, a character string can also be provided, specifying a FASTA file containing all the reference sequences in a genome.
The latter may be necessary to synchronise the fragments with the genome used for alignment.
}

\author{Aaron Lun}

\examples{
require(BSgenome.Ecoli.NCBI.20080805)

cutGenome(Ecoli, "AAGCTT", overhang=4L) # HindIII
cutGenome(Ecoli, "CCGCGG", overhang=2L) # SacII
cutGenome(Ecoli, "AGCT", overhang=0L) # AluI

# Trying with FastA files.
x <- system.file("extdata", "fastaEx.fa", package="Biostrings")
cutGenome(x, "AGCT", overhang=2)
cutGenome(x, "AGCT", overhang=4)
}

\seealso{
\code{\link{matchPattern}}
}

