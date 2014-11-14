#' @export allelic.richness
#' 
#' @title Allelic Richness 
#' @description Calculate allelic richness for diploid data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param use.na.rows logical. FALSE deletes any row with an NA for either allele. 
#'   TRUE uses all unique (non-NA) alleles in frequency.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

allelic.richness <- function(g, use.na.rows = FALSE) {  
  stopifnot.gtypes(g, "diploid")
  
  sapply(attr(g, "locus.names"), function(x) {
    locus <- g$genotypes[, locus.cols(x, g)]
    locus[locus == -1] <- NA
    locus <- if(use.na.rows) as.vector(na.omit(locus)) else na.omit(as.vector(locus))
    length(unique(locus)) / (length(locus) / 2)
  })
}