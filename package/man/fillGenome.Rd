\name{fillGenome}
\alias{fillGenome}

\title{Fill-in genomic restriction sites}

\description{Perform restriction site cleavage, blunt-end repair and ligation for a target genome.}

\usage{
fillGenome(bs, pattern, outfile, overhang=4L, spacer=50L)
}

\arguments{
	\item{bs}{a \code{BSgenome} object}
	\item{pattern}{character string describing the recognition site}
	\item{outfile}{character string specifying the name of the output FASTA file}
	\item{overhang}{integer scalar specifying the length of the 5' overhang}
	\item{spacer}{integer scalar specifying the length of a spacer sequence to insert between blunt ends}
}

\value{
A \code{GRangesList} object containing elements \code{original} and \code{modified}. The former describes the 
restriction fragment boundaries on the original genome, the latter describes the same boundaries on the modified genome.
}

\details{
This function simulates a restriction digestion using \code{\link{cutGenome}}.
It performs blunt-end ligation by filling in the 5' overhangs. The blunt ends
are ligated back together with an intervening spacer sequence. This is designed
to mimic the manipulations in the Hi-C procedure. The spacer discourages the
incorrect alignment of chimeric reads as non-chimeric reads spanning
restriction sites.

The sequences of the modified chromosomes are saved in FASTA format in the
specified \code{outfile}. This can be used to construct an index for alignment
of Hi-C reads. The boundaries of the restriction fragments on each chromosome
are also provided as output. This can be saved to file and used to direct
alignment of chimeric reads e.g. splitting soft-clipped reads according to
known breakpoints. 
}

\author{Aaron Lun}

\examples{
require(BSgenome.Ecoli.NCBI.20080805)
fillGenome(Ecoli, "AAGCTT", "out.fa", overhang=4L)
fillGenome(Ecoli, "AAGCTT", "out.fa", overhang=2L)
fillGenome(Ecoli, "AAGCTT", "out.fa", spacer=10L)
}

\seealso{
\code{\link{cutGenome}}
}

