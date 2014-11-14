#' @export as.data.frame.gtypes
#' 
#' @title Coerce to a Data Frame
#' @description Coerce a gtypes object to a data.frame.
#'
#' @usage \method{as.data.frame}{gtypes}(x, ...)
#' 
#' @param x a \code{\link{gtypes}} object.
#' @param ... additional arguments to be passed to or from methods.
#' 
#' @return a data.frame with the columns 'id', 'strata', and two columns for every 
#'   diploid locus, or one column of haplotypes.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' @seealso \code{link{as.matrix.gtypes}}

as.data.frame.gtypes <- function(x, ...) as.data.frame(as.matrix(x), stringsAsFactors = FALSE)