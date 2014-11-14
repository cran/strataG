#' @export stat.chi2
#' 
#' @title CHI-squared Analysis of Population Structure
#' 
#' @param g a \code{gtypes} object.
#' 
#' @return a \code{\link{gtype.struct.stat}} list with the CHI-squared estimate.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.chi2(msats)

stat.chi2 <- function(g) {
  g <- check.gtypes(g)
  
  chi2.func <- if(is.diploid(g)) chi2D_C else chi2H_C
  est <- chi2.func(g$genotypes[, -1, drop = FALSE], g$genotypes[, "strata"])
  
  result <- list(stat.name = "Chi2", estimate = est, strata.freq = table(decode.strata(g)))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.chi2) <- c(class(stat.chi2), "gtype.struct.func")