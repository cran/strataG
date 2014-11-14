#' @export strata.split
#' 
#' @title Split Strata
#' @description Return a list of gtypes for each strata
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param strata a character vector giving a subset of strata to select. If \code{NULL} then a list with 
#'   all strata is created.
#' @param keep.unstratified logical. Create a gtypes object with unstratified individuals? 
#' @param ... other parameters to be passed to \code{\link{subset.gtypes}}
#' 
#' @return A named list where each element is a \code{gtypes} object for a single stratum in \code{g}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{subset.gtypes}} to create a single \code{gtypes} object for specific strata.

strata.split <- function(g, strata = NULL, keep.unstratified = FALSE, ...) {
  stopifnot.gtypes(g)
  
  if(!is.null(strata)) g <- subset(g, strata = strata)
  
  strata.vec <- attr(g, "strata.names")
  if(!keep.unstratified) strata.vec <- na.omit(strata.vec)
  result <- sapply(strata.vec, function(st) subset(g, strata = st, ...), simplify = FALSE)
  names(result) <- strata.vec
  result
}