#' @export stat.gst.dbl.prime
#' 
#' @title G''st Analysis of Population Structure
#' 
#' @param g a \code{gtypes} object
#' 
#' @details Calculate G''st from Meirmans and Hedrick 2010
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.gst.dbl.prime(msats)

stat.gst.dbl.prime <- function(g) {
  g <- check.gtypes(g, "diploid")
  
  hets <- hStats_C(g$genotypes[, -1, drop = FALSE], g$genotypes[, "strata"])
  est <- gstDblPrime_C(g$genotypes[, "strata"], hets)
  if(is.nan(est)) est <- NA
  
  result <- list(stat.name = "G''st", estimate = est, strata.freq = table(decode.strata(g)))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.gst.dbl.prime) <- c(class(stat.gst.dbl.prime), "gtype.struct.func")