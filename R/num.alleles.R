#' @name num.alleles
#' @export num.alleles
#' 
#' @title Number of Alleles
#' @description Return the number of alleles for each locus.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param use.na.rows logical. FALSE deletes any row with an NA for either allele. 
#'   TRUE uses all unique (non-NA) alleles in frequency.
#' 
#' @return vector of number of alleles per locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

num.alleles <- function(g, use.na.rows = TRUE) {  
  stopifnot.gtypes(g)
  
  if(is.haploid(g)) {
    length(unique(decode.loci(g)))
  } else {
    num.alleles <- sapply(1:length(attr(g, "locus.names")), function(i) {
      locus <- g$genotypes[, locus.cols(i)]
      locus[locus == -1] <- NA
      locus <- if(use.na.rows) as.vector(na.omit(locus)) else na.omit(as.vector(locus))
      length(unique(locus))
    })
    names(num.alleles) <- attr(g, "locus.names")
    num.alleles
  }
}