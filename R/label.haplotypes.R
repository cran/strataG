#' @rdname label.haplotypes
#' @export label.haplotypes label.haplotypes.gtypes label.haplotypes.default
#' 
#' @title Find and label haplotypes
#' @description Identify and group sequences that share the same haplotype.
#' 
#' @param x list of DNA sequences or a \code{\link{gtypes}} object.
#' @param prefix a character string giving prefix to be applied to numbered haplotypes.
#'   If NULL, haplotypes will be labeled with first label from original sequences.
#' @param ignore.gaps logical. Remove sites with gaps when comparing sequences?
#' @param ... Further arguments to be passed to \code{label.haplotypes.default}.
#' 
#' @note If any sequences contain missing data (N's), a more rigorous matching algorithm is 
#'   conducted and can take longer, especially if there are many samples or long sequences.
#' 
#' @return List of haplotypes and assignmnents or gtypes object with new haplotypes.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

label.haplotypes <- function(x, prefix = NULL, ignore.gaps = FALSE) UseMethod("label.haplotypes")