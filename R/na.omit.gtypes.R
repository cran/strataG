#' @name na.omit.gtypes
#' @export na.omit.gtypes
#' 
#' @title Omit Missing Data
#' @description Omit samples or loci with missing data in \code{gtypes} object.
#' 
#' @usage \method{na.omit}{gtypes}(object, strata = TRUE, loci = TRUE, ...)
#' 
#' @param object a \code{\link{gtypes}} object.
#' @param strata logical. Remove unstratified samples?
#' @param loci logical. Remove samples with missing data for all loci?
#' @param ... additional arguments to be passed to or from methods. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

na.omit.gtypes <- function(object, strata = TRUE, loci = TRUE, ...) {
  stopifnot.gtypes(object)
  
  g <- decode(object)
  missing.strata <- if(strata) which(is.na(g$genotypes[, "strata"])) else NULL
  missing.loci <- if(loci) which(apply(g$genotypes[, -1, drop = FALSE], 1, function(x) all(is.na(x)))) else NULL
  to.remove <- unique(c(missing.strata, missing.loci))
  to.keep <- rownames(g$genotypes)[-to.remove]
  if(length(to.remove) > 0) subset(object, ids = to.keep) else object
}