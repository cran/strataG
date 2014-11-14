#' @export as.matrix.gtypes
#' 
#' @title Coerce to a Matrix
#' @description Coerce a gtypes object to a matrix.
#' 
#' @usage \method{as.matrix}{gtypes}(x, ...)
#' 
#' @param x a \code{\link{gtypes}} object.
#' @param ... additional arguments to be passed to or from methods.
#' 
#' @return A matrix with the columns 'id', 'strata', and two columns for every 
#'   diploid locus, or one column of haplotypes
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{as.data.frame.gtypes}}

as.matrix.gtypes <- function(x, ...) {
  stopifnot.gtypes(x)
  
  gtypes.mat <- cbind(rownames(x$genotypes), decode.strata(x), decode.loci(x))
  locus.colnames <- if(is.diploid(x)) {
    paste(rep(attr(x, "locus.names"), each = 2), 1:2, sep = ".")
  } else {
    attr(x, "locus.names")
  }
  colnames(gtypes.mat) <- c("id", "strata", locus.colnames)
  rownames(gtypes.mat) <- NULL
  gtypes.mat
}