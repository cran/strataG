#' @export hap.freqs
#' 
#' @title Haplotype Frequencies
#' @description Create a table of haplotype frequencies in each strata.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @return A table of haplotype frequencies.
#' 
#' @note Function assumes that haplotypes have been already labelled. If not, use 
#'   \code{\link{label.haplotypes}} first.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.haps)
#' 
#' mtdna <- gtypes(dolph.strata, id.col = 1, strata.col = 2, locus.col = 4, dna.seq = dolph.haps)
#' hap.freqs(mtdna)

hap.freqs <- function(g) {
  stopifnot.gtypes(g, "haploid")
  
  g <- decode(g)
  table(Haplotype = g$genotypes[, 2], Strata = g$genotypes[, "strata"])
}