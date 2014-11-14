#' @export stat.fst
#' 
#' @title Fst Analysis of Population Structure
#' 
#' @param g a \code{gtypes} object.
#' 
#' @return a \code{\link{gtype.struct.stat}} list with the Fst estimate.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Weir, B.S., Cockerham, C.C. 1984. Estimating F-statistics
#'  for the analysis of population structure. Evolution 38(6):1358-1370.
#'  
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.fst(msats)

stat.fst <- function(g) {  
  g <- check.gtypes(g)
  
  if(is.haploid(g)) {
    g$sequences <- NULL
    return(stat.phist(g))
  }

  est <- fst_C(g$genotypes[, -1, drop = FALSE], g$genotypes[, "strata"])
  if(is.nan(est)) est <- NA
  
  result <- list(stat.name = "Fst", estimate = est, strata.freq = table(decode.strata(g)))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.fst) <- c(class(stat.fst), "gtype.struct.func")