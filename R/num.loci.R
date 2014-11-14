#' @rdname locus.names

num.loci <- function(locus.data) { 
  if(!(is.matrix(locus.data) | is.data.frame(locus.data))) {
    stop("'locus.data' is not a matrix or data.frame.")
  }
  if(ncol(locus.data) == 1) return(1)
  if(ncol(locus.data) %% 2) {
    stop("'locus.data' has more than one column and an odd number of columns.")
  }
  ncol(locus.data) / 2
}