#' @export variable.sites
#' 
#' @title Variable Sites
#' @description Identify variable sites among sequences.
#'  
#' @param x a \code{\link{gtypes}} object with sequences, a list of sequences, or a consensus sequence.
#'   Sequences must be aligned. 
#' @param bases character vector of bases to consider.
#' 
#' @return A list with: \tabular{ll}{
#'   \code{site} \tab a list of sequences composed of variable sites. \cr
#'   \code{site.freqs} \tab a matrix of base pair frequencies by site. \cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{fixed.sites}}

variable.sites <- function(x, bases = c("a", "c", "g", "t", "-")) {
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x)
  stopifnot.aligned(x)
  
  site.freqs <- base.freqs(x, bases)$site.freqs
  var.site <- apply(site.freqs, 2, function(site) sum(site > 0) > 1)
  var.seqs <- lapply(x, function(x) tolower(x[var.site]))
  var.site.freqs <- as.matrix(site.freqs[, var.site], nrow = nrow(site.freqs))
  colnames(var.site.freqs) <- which(var.site)
  list(sites = var.seqs, site.freqs = var.site.freqs)
}