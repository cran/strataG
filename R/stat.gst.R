#' @export stat.gst
#' 
#' @title Gst Analysis of Population Structure
#' 
#' @param g a \code{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Nei 1973
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.gst(msats)


stat.gst <- function(g) {
  g <- check.gtypes(g, "diploid")
  
  hets <- hStats_C(g$genotypes[, -1, drop = FALSE], g$genotypes[, "strata"])
  est <- gst_C(g$genotypes[, "strata"], hets)
  if(is.nan(est)) est <- NA
  
  result <- list(stat.name = "Gst", estimate = est, strata.freq = table(decode.strata(g)))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.gst) <- c(class(stat.gst), "gtype.struct.func")