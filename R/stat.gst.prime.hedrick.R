#' @export stat.gst.prime.hedrick
#' 
#' @title Hedrick's G'st Analysis of Population Structure
#' 
#' @param g a \code{gtypes} object.
#' 
#' @details Calculate G'st from Hedrick 2005
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.gst.prime.hedrick(msats)

stat.gst.prime.hedrick <- function(g) {
  g <- check.gtypes(g, "diploid")
  
  hets <- hStats_C(g$genotypes[, -1, drop = FALSE], g$genotypes[, "strata"])
  gst <- gst_C(g$genotypes[, "strata"], hets)
  est <- gstPrimeHedrick_C(g$genotypes[, "strata"], hets, gst)
  if(is.nan(est)) est <- NA

  result <- list(stat.name = "G'st (Hedrick 2005)", estimate = est, strata.freq = table(decode.strata(g)))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.gst.prime.hedrick) <- c(class(stat.gst.prime.hedrick), "gtype.struct.func")