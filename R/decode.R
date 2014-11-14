#' @rdname decode
#' @export decode decode.strata decode.loci decode.sequences
#' 
#' @title Decode \code{gtypes}
#' @description Functions to decode \code{gtypes} objects.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @return \tabular{ll}{
#'   \code{decode} \tab a list with the results from \code{decode.strata}, \code{decode.locus.data}, and \code{decode.sequences}.\cr
#'   \code{decode.strata} \tab a named vector of strata designations.\cr
#'   \code{decode.locus.data} \tab a matrix of haplotypes or genotypes.\cr
#'   \code{decode.sequences} \tab a list of DNA sequences if present, otherwise \code{NULL}.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

decode <- function(g) {
  result <- list(genotypes = cbind(strata = decode.strata(g), decode.loci(g)),
       sequences = decode.sequences(g)
  )
  attr(result, "locus.names") <- attr(g, "locus.names")
  attr(result, "description") <- attr(g, "description")
  result
}