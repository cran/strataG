#' @export create.consensus
#' 
#' @title Consensus Sequence
#' @description Return a consensus sequence from set of aligned sequences, introducing IUPAC ambiguity
#'   codes where necessary.
#' 
#' @param x a \code{\link{gtypes}} object with aligned sequences, or a list of aligned DNA sequences.
#' @param ignore.gaps logical. Ignore gaps at a site when creating consensus. If true, then bases with a 
#'   gap are removed before consensus is calculated. If false and a gap is present, then the result is a gap.
#' 
#' @return A character vector of the consensus sequence.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.seqs)
#' create.consensus(dolph.seqs)

create.consensus <- function(x, ignore.gaps = FALSE) { 
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x)
  stopifnot.aligned(x)
  x <- do.call(rbind, lapply(x, tolower))
  apply(x, 2, iupac.code, ignore.gaps = ignore.gaps)
}