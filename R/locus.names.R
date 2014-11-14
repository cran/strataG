#' @rdname locus.names
#' @export locus.names locus.cols num.loci
#' 
#' @title Locus Information
#' @description Helper functions for locus information from \code{gtypes} objects.
#' 
#' @param locus.data vector or matrix of genetic data.
#' @param x vector of locus numbers or locus names.
#' @param g a \code{\link{gtypes}} object to look for locus names in if \code{x} is a character.
#' 
#' @return \tabular{ll}{
#'   \code{locus.names} \tab a character vector of names of loci.\cr
#'   \code{locus.cols} \tab a matrix of columns in the \code{genotypes} element of \code{g} for the specified loci in \code{x}.\cr
#'   \code{num.loci} \tab x.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

locus.names <- function(locus.data) {  
  if(is.vector(locus.data)) {
    locus.data <- cbind(locus.data)
    colnames(locus.data) <- NULL
  }
  n.loc <- num.loci(locus.data)
  if(n.loc == 1) return("haplotype")
  
  if(is.null(colnames(locus.data))) {
    # return generic names if no colnames assigned
    locus.nums <- formatC(1:n.loc, digits = floor(log10(n.loc)), flag = "0")  
    paste("Locus_", locus.nums, sep = "")
  } else {
    sapply(seq(1, ncol(locus.data), by = 2), function(col1) {
      loc1 <- unlist(strsplit(colnames(locus.data)[col1], ""))
      loc2 <- unlist(strsplit(colnames(locus.data)[col1 + 1], ""))
      same.char <- loc1 == loc2
      # return first allele colname if all characters are different or all the same
      if(all(same.char) || all(!same.char)) return(colnames(locus.data)[col1])
      # otherwise return only equivalent portion of colnames minus trailing non-alphanumeric characters
      loc.name <- loc1[same.char]
      while(!(loc.name[length(loc.name)] %in% c(LETTERS, letters, 0:9))) loc.name <- loc.name[-length(loc.name)]
      paste(loc.name, collapse = "")
    })
  }
}
