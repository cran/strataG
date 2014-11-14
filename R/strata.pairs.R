#' @export strata.pairs
#' 
#' @title Strata Pairs
#' @description Make a data.frame of all pairs of strata.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

strata.pairs <- function(g) {
  strata <- sort(unique(decode.strata(g)))
  if(length(strata) < 2) stop("'g' has fewer than 2 strata")
  df <- as.data.frame(t(combn(strata, 2)), stringsAsFactors = FALSE)
  colnames(df) <- c("strata.1", "strata.2")
  df
}