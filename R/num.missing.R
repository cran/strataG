#' @export num.missing
#' 
#' @title Number Missing Data
#' @description Calculate the number of individuals with missing data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param prop logical determining whether to return proportion missing.
#' 
#' @return a vector of loci with number (or, if \code{prop = TRUE}, the proportion) 
#'   of individuals missing data for at least one allele.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' num.missing(msats)
#' num.missing(msats, prop = TRUE)

num.missing <- function(g, prop = FALSE) {  
  stopifnot.gtypes(g)
  
  num.missing <- if(is.haploid(g)) {
    sum(g$genotypes[, 2] == -1)
  } else {
    sapply(attr(g, "locus.names"), function(x) {
      sum(apply(g$genotypes[, locus.cols(x, g), drop = FALSE], 1, function(x) any(x == -1)))
    })
  }
  names(num.missing) <- attr(g, "locus.names")
  if(prop) num.missing / nrow(g$genotypes) else num.missing
}