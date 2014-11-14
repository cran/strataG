#' @rdname is.diploid
#' @export is.diploid is.haploid
#' 
#' @title Check Ploidy
#' @description Check if a gtypes object contains haploid or diploid data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @details Checks if the \code{locus.data} element of \code{g} is one column (haploid) 
#'   or an even-numbered column (diploid) matrix.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

is.diploid <- function(g) {
  stopifnot.gtypes(g)
  ncol(g$genotypes) > 2
}