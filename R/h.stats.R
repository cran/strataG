#' @export h.stats
#' 
#' @title Estimate Fixation Indices
#' @description Estimate Ho, Hs, and Ht from diploid data.
#' 
#' @param g a \code{gtypes} object.
#' 
#' @return matrix of Ho, Hs, and Ht for each locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Nei, M. and R.K. Chesser. 1983. Estimation of fixation indices and gene diversities. 
#'   Ann. Hum. Genet. 47:253-259.

h.stats <- function(g) {  
  stopifnot.gtypes(g, "diploid")
  
  result <- hStats_C(g$genotypes[, -1], g$genotypes[, "strata"])
  h.labels <- c("Ho", "Hs", "Ht")
  if(is.matrix(result)) rownames(result) <- h.labels else names(result) <- h.labels
  result
}