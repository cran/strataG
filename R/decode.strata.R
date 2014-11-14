#' @rdname decode

decode.strata <- function(g) {
  stopifnot.gtypes(g)
  decoded <- attr(g, "strata.names")[as.character(g$genotypes[, "strata"])]
  names(decoded) <- rownames(g$genotypes)
  decoded
}