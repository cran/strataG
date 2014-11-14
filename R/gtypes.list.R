#' @rdname gtypes

gtypes.list <- function(gen.data, strata.vec, ... ) {
  if(!is.character(strata.vec) & !is.numeric(strata.vec) & length(strata.vec) != length(gen.data)) {
    stop("'strata.vec' must be a character or numeric vector the same length as gen.data")
  }
  gen.mat <- cbind(id = names(gen.data), strata = strata.vec, haplotypes = names(gen.data))
  gtypes.default(gen.mat, id.col = 1, strata.col = 2, locus.col = 3, dna.seq = gen.data, ...) 
}
