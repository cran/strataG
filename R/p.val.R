#' @export p.val
#' 
#' @title Permutation Test P-value
#' @description Calculate the p-value for a permutation test.
#' 
#' @param obs observed value.
#' @param null.dist vector of values from permutation null distribution.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

p.val <- function(obs, null.dist) (sum(null.dist >= obs) + 1) / (length(null.dist) + 1)