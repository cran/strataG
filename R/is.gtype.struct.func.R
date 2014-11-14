#' @title Check Population Structure Function
#' @description Check if function is a valid population structure function. Used
#'   in overall and pairwise test functions.
#' 
#' @param x an \code{R} object.
#' @details checks if \code{x} is a \code{\link{gtypes}} population structure function. 
#'   Used in permutation test functions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

is.gtype.struct.func <- function(x) is.function(x) & "gtype.struct.func" %in% class(x)
