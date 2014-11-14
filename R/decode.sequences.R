#' @rdname decode

decode.sequences <- function(g) {
  stopifnot.gtypes(g)
  
  if(!is.null(g$sequences)) {
    decoded <- g$sequences
    names(decoded) <- attr(g, "allele.names")[[1]][names(g$sequences)]
    decoded
  } else NULL
}