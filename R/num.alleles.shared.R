#' @name num.alleles.shared
#' @export num.alleles.shared
#' 
#' @title Number of Shared Alleles
#' @description Calculate the number of alleles shared between two samples at each locus. 
#' 
#' @param id1,id2 ids of samples to compare.
#' @param g a \code{\link{gtypes}} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

num.alleles.shared <- function(id1, id2, g) {
  sapply(attr(g, "locus.names"), function(locus) {
    lc <- locus.cols(locus, g)
    gtype1 <- g$genotypes[id1, lc]    
    gtype2 <- g$genotypes[id2, lc]
    if(any(c(gtype1, gtype2) == -1)) return(NA)
    sum(gtype1 %in% gtype2) + sum(gtype2 %in% gtype1) 
  })
}