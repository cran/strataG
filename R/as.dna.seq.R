#' @rdname as.dna.seq
#' @export as.dna.seq is.aligned.seq is.dna.seq stopifnot.aligned
#' 
#' @title Format List of DNA Sequences
#' @description Format and test DNA sequences stored in lists.
#' 
#' @param x an R object of either matrix, data.frame, character vector, list, or DNAbin.
#' @param show.warnings logical - show warnings for is.dna.seq describing check failures?
#' 
#' @details \tabular{ll}{
#'   \code{as.dna.seq} \tab converts a \code{matrix}, \code{data.frame}, \code{vector}, 
#'     or \code{\link[ape]{DNAbin}} object to a list of character vectors. Nucleotides are 
#'     stored as lower-case.\cr
#'   \code{is.aligned.seq} \tab tests if object is an aligned set of DNA sequences. 
#'     Each sequence must be the same length.\cr
#'   \code{is.dna.seq} \tab tests if object is a list where each element is a character 
#'     vector and every element in each vector is a valid IUPAC nucleotide code.\cr
#'   \code{stopifnot.aligned} \tab stops execution if object is not an aligned set of DNA sequences.\cr
#' }
#'  
#' @author Eric Archer \email{eric.archer@@noaa.gov}
 
as.dna.seq <- function(x) {
  x.names <- NULL
  if(any(c("matrix", "data.frame") %in% class(x))) {
    if(is.null(rownames(x))) rownames(x) <- 1:nrow(x)
    x.names <- rownames(x)
    x <- sapply(rownames(x), function(i) as.vector(x[i, ]), simplify = FALSE)
  } else if(any(c("character", "DNAbin")  %in% class(x))) {
    if(is.null(names(x))) names(x) <- 1:nrow(x)
    x.names <- names(x)
    x <- strsplit(as.character(x), "")
  } else if("list" %in% class(x)) x.names <- names(x)
  
  if("list" %in% class(x)) {
    x <- lapply(x, function(y) tolower(as.character(y)))
    class(x) <- c("list", "dna.seq")
    names(x) <- if(is.null(x.names)) 1:length(x) else x.names
    return(x)
  } else stop("'x' is of unknown type")
}

