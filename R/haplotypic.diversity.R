#' @export haplotypic.diversity
#' 
#' @title Haplotypic Diversity
#' @description Calculate haplotypic diversity of all samples.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @note Function assumes that haplotypes in \code{g} have been already labelled. If not, use \code{\link{label.haplotypes}} first.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}  
#' 
#' @seealso \code{\link{diversity}}, \code{\link{nucleotide.diversity}}

haplotypic.diversity <- function(g) {
  stopifnot.gtypes(g, "haploid") 
  
  haps <- g$genotypes[, 2]
  haps[haps == -1] <- NA
  diversity(haps)  
}