#' @export obsvd.het 
#' 
#' @title Observed Heterozygosity 
#' @description Calculate observed heterozygosity for diploid data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{exptd.het}}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' obsvd.het(msats)

obsvd.het <- function(g) {
  stopifnot.gtypes(g, "diploid")
  
  sapply(attr(g, "locus.names"), function(x) {
    locus <- g$genotypes[, locus.cols(x, g), drop = FALSE]
    locus[locus == -1] <- NA
    locus <- na.omit(locus)
    heterozygote <- ifelse(locus[, 1] != locus[, 2], TRUE, FALSE)
    sum(heterozygote) / nrow(locus)
  })
}