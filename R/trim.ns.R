#' @export trim.ns
#' 
#' @title Trim N's From Sequences
#' @description Removes N's from beginning and end of sequences.
#' 
#' @param dna.seq a list of DNA sequences.
#' 
#' @return a list of DNA sequences.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

trim.ns <- function(dna.seq) {
  stopifnot.aligned(dna.seq)
  lapply(dna.seq, function(x) {
    x.vec <- if(is.vector(x)) x else x[1, ]
    x.vec <- paste(x.vec, collapse = "")
    start <- gregexpr("^[n]+", x.vec)[[1]]
    end <- gregexpr("[n]+$", x.vec)[[1]]
    start <- ifelse(start == -1, 1, attr(start, "match.length") + 1)
    end <- ifelse(end == -1, nchar(x.vec), end)
    if(is.vector(x)) x[start:end] else x[, start:end]
  })
}