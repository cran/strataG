#' @rdname as.dna.seq

stopifnot.aligned <- function(x, show.warnings = TRUE) {
  if(!is.aligned.seq(x, show.warnings = show.warnings)) stop("'x' is not aligned set of sequences.")
}