#' @export theta
#' @importFrom pegas theta.h
#' 
#' @title Theta
#' @description Calculate theta from heterozygosity of each locus.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @return vector of theta values for each locus.
#' 
#' @details Calculates theta for each locus using the \code{\link[pegas]{theta.h}} function.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' theta(msats)

theta <- function(g) {
  stopifnot.gtypes(g)
  
  theta.locus <- if(is.haploid(g)) {
    theta.h(as.factor(na.omit(g$genotypes[, 2])))
  } else {
    sapply(1:length(attr(g, "locus.names")), function(i) {
      locus <- na.omit(g$genotypes[, locus.cols(i)])
      theta.h(as.factor(as.vector(locus)))
    })
  }
  names(theta.locus) <- attr(g, "locus.names")
  theta.locus
}