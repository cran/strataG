#' @export restratify
#' 
#' @title Restratify gtypes
#' @description Change the stratification of samples in a \code{\link{gtypes}} object.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param new.strata a character or numeric vector specifying new strata for samples. 
#'   If vector is named, names are assumed to be sample IDs. If not named, vector must be same length as number
#'   of samples in \code{g}.
#'   
#' @return new \code{\link{gtypes}} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

restratify <- function(g, new.strata) {
  stopifnot.gtypes(g)
  
  # check format of 'new.strata'
  if(!is.vector(new.strata)) stop("'new.strata' must be a character or numeric vector")
  new.strata <- as.character(new.strata[!is.na(new.strata)])
  if(length(new.strata) != nrow(g$genotypes)) {
    if(is.null(names(new.strata))) stop("'new.strata' must be named with sample IDs if it does not have the same number of samples as 'g'")
  } else if(is.null(names(new.strata))) names(new.strata) <- rownames(g$genotypes)
  
  # check if IDs in 'new.strata' are in 'g'
  not.found <- setdiff(names(new.strata), rownames(g$genotypes))
  if(length(not.found) == length(new.strata)) stop("No ids in 'g' gtypes found in 'new.strata' vector")
  if(length(not.found) > 0) {
    not.found <- paste(not.found, collapse = ", ")
    warning(paste("The following ids in 'new.strata' vector are not in 'g' gtypes:", not.found))
  }
  
  # replace strata of IDs found with 'new.strata'
  found <- intersect(names(new.strata), rownames(g$genotypes))
  new.strata <- new.strata[found]
  desc <- attr(g, "description")
  g <- decode(g)
  gen.mat <- cbind(id = rownames(g$genotypes), g$genotypes)
  gen.mat[found, "strata"] <- new.strata[found]
  gtypes(gen.mat, dna.seq = g$sequences, description = desc)
}