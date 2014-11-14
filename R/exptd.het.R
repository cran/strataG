#' @export exptd.het
#' 
#' @title Expected Heterozygosity
#' @description Calculate expected heterozygosity for diploid data.
#' 
#' @param g a \code{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{obsvd.het}}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' exptd.het(msats)

exptd.het <- function(g) {  
  stopifnot.gtypes(g, "diploid")
  
  sapply(attr(g, "locus.names"), function(x) {
    locus <- as.vector(g$genotypes[, locus.cols(x, g)])
    locus[locus == -1] <- NA
    diversity(locus)
  })
}
