#' @rdname is.diploid

is.haploid <- function(g) {
  stopifnot.gtypes(g)
  ncol(g$genotypes) == 2
}