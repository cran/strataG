#' @export pct.unique.haplotypes
#' 
#' @title Percent Unique Haplotypes
#' @description Calculate the percent of haplotypes that are unique.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

pct.unique.haplotypes <- function(g) {  
  stopifnot.gtypes(g, "haploid")
  
  haps <- g$genotypes[, 2]
  haps[haps == -1] <- NA
  sum(table(haps) == 1) / nrow(g$genotypes)
}