#' @title Check gtype Population Structure Function
#' @description Stop execution if function is not a valid population structure
#' 
#' @param func function to test.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

stopifnot.gtype.struct.func <- function(func) if(!is.gtype.struct.func(func)) stop("Improper 'gtype.struct.func' supplied")
