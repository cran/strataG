#' @rdname as.dna.seq

is.aligned.seq <- function(x, show.warnings = FALSE) {
  if(!is.dna.seq(x, show.warnings = show.warnings)) return(FALSE)
  if(length(unique(sapply(x, length))) != 1) {
    if(show.warnings) warning("All sequences in 'x' are not the same length.")
    return(FALSE)
  }
  TRUE
}