#' @export fixed.sites
#' 
#' @title Fixed Sites
#' @description Identify fixed sites among sequences.
#' 
#' @param x a \code{\link{gtypes}} object with sequences, a list of sequences, or a consensus sequence. 
#'   Sequences must be aligned.
#' @param bases a character vector of valid bases to consider.
#' 
#' @author Eric Archer <eric.archer@@noaa.gov>
#' 
#' @seealso \code{\link{variable.sites}}

fixed.sites <- function(x, bases = c("a", "c", "g", "t", "-")) {
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x)
  stopifnot.aligned(x)
  
  seq.mat <- tolower(do.call(rbind, x))
  bases <- tolower(bases)
  is.fixed <- sapply(1:ncol(seq.mat), function(i) {
    bp <- seq.mat[, i]
    bp <- bp[bp %in% bases]
    length(unique(bp)) == 1
  })
  sites <- which(is.fixed)
  fixed.bp <- sapply(sites, function(i) {
    bps <- seq.mat[, i]
    unique(bps[bps %in% bases])[1]
  })
  names(fixed.bp) <- sites
  fixed.bp
}