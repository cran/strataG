#' @export base.freqs
#' 
#' @title Base Frequencies
#' @description Calculate nucleotide base frequencies along a sequence.
#' 
#' @param x a \code{\link{gtypes}} object with aligned sequences or a list of aligned DNA sequences.
#' @param bases character vector of bases. Must contain valid IUPAC codes. If NULL, will return summary 
#'   of frequencies of observed bases.
#' 
#' @return list containing:
#' \tabular{ll}{
#'   \code{site.freqs} \tab a matrix of base frequencies at each site.\cr
#'   \code{base.freqs} \tab a vector of overall base frequencies.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

base.freqs <- function(x, bases = NULL) {
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x$sequences)
  stopifnot.aligned(x)
  
  seq.mat <- do.call(rbind, lapply(x, tolower))
  bases <- if(is.null(bases)) {
    rownames(iupac.mat)
  } else {
    tolower(as.character(bases))
  }
  site.freqs <- apply(seq.mat, 2, function(site) table(factor(site, levels = bases)))
  colnames(site.freqs) <- 1:ncol(site.freqs)
  base.freqs <- table(factor(as.vector(seq.mat), levels = bases))
  base.freqs <- base.freqs / sum(base.freqs)
  
  list(site.freqs = site.freqs, base.freqs = base.freqs)
}