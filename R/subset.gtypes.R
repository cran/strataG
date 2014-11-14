#' @export subset.gtypes
#' 
#' @title Subset a gtypes Object
#' @description Return a subset of specific individuals, strata, or loci from a gtypes object.
#' 
#' @usage \method{subset}{gtypes}(x, ids, strata, loci, remove.sequences, ...)
#' 
#' @param x a \code{\link{gtypes}} object.
#' @param ids a character vector specifying which ids to keep.
#' @param strata a character vector specifying which strata to keep.
#' @param loci a character vector specifying which loci to keep.
#' @param remove.sequences logical. For haploid objects, remove sequences for samples that are excluded?
#' @param ... other parameters (ignored). 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

subset.gtypes <- function(x, ids = NULL, strata = NULL, loci = NULL, remove.sequences = FALSE, ...) {
  stopifnot.gtypes(x)
  g <- decode(x)
  
  if(!is.null(ids)) {
    if(!(is.atomic(ids) & is.character(ids))) stop("'ids' must be a character vector.")
    ids <- unique(ids)
    ids.found <- ids %in% rownames(g$genotypes)
    if(!all(ids.found)) {
      missing.ids <- paste(rownames(g$genotypes)[!(ids.found)], collapse = ", ")
      stop("The following ids could not be found:", missing.ids)
    }
  } else ids <- rownames(g$genotypes)
  
  if(!is.null(strata)) {
    if(!(is.atomic(strata) & is.character(strata))) stop("'strata' must be a character vector.")
    strata <- unique(strata)
    strata.found <- strata %in% g$genotypes[, "strata"]
    if(!all(strata.found)) {
      missing.strata <- paste(strata[!(strata.found)], collapse = ", ")
      stop("The following strata could not be found:", missing.strata)
    }
    ids <- intersect(ids, rownames(g$genotypes)[g$genotypes[, "strata"] %in% strata])
  }
  
  locus.names <- attr(g, "locus.names")
  locus.nums <- if(!is.null(loci) & ncol(g$genotypes) > 2) {
    if(!(is.atomic(loci) & is.character(loci))) stop("'loci' must be a character vector.")
    loci.found <- loci %in% locus.names
    if(!all(loci.found)) {
      missing.loci <- paste(loci[!(loci.found)], collapse = ", ")
      stop(paste("The following loci could not be found:", missing.loci))
    }
    which(locus.names %in% loci)
  } else 1:length(locus.names)
  loc.cols <- if(ncol(g$genotypes) == 2) 2 else as.vector(locus.cols(locus.nums))
    
  gen.mat <- cbind(ids = ids, g$genotypes[ids, c(1, loc.cols), drop = FALSE])
  dna.seq <- g$sequences
  if(!is.null(g$sequences) & remove.sequences) dna.seq <- g$sequences[unique(g$genotypes[, 2])]

  gtypes(gen.mat, dna.seq = dna.seq, description = attr(x, "description"))
}