% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lowFreqSubs.R
\name{lowFreqSubs}
\alias{lowFreqSubs}
\title{Low Frequency Substitutions}
\usage{
lowFreqSubs(x, min.freq = 3, motif.length = 10, simplify = TRUE)
}
\arguments{
\item{x}{a \code{\link[ape]{DNAbin}} object.}

\item{min.freq}{minimum frequency of base to be flagged.}

\item{motif.length}{length of motif around low frequency base to output.}

\item{simplify}{if there is a single locus, return result in a simplified
form? If \code{FALSE} a list will be returned wth one element per locus.}
}
\value{
data.frame listing id, site number, and motif around low frequency 
  base call.
}
\description{
Check nucleotide sites for low frequency substitutions.
}
\examples{
data(dolph.haps)

lowFreqSubs(dolph.haps)

}
\author{
Eric Archer \email{eric.archer@noaa.gov}
}
