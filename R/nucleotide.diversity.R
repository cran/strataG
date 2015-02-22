#' @rdname nucleotide.diversity
#' @export nucleotide.diversity
#'   
#' @title Nucleotide Diversity
#' @description Calculate nucleotide diversity for set of haplotypes.
#' 
#' @param x a list of sequences or a \code{\link{gtypes}} object.
#' @param bases nucleotides to consider when calculating diversity.
#' 
#' @return Nucleotide diversity by site. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

nucleotide.diversity <- function(x, bases = c("a", "c", "g", "t")) {
  seq.mat <- if(is.gtypes(x, show.warnings = FALSE)) {
    stopifnot.gtypes(x, "haploid")
    stopifnot.aligned(x$sequences)
    hap.freq <- table(factor(x$genotypes[, 2], levels = names(x$sequences)))
    seq.mat <- tolower(do.call(rbind, x$sequences))
    seq.mat[rep(1:nrow(seq.mat), hap.freq), ]
  } else {
    stopifnot.aligned(x)
    tolower(do.call(rbind, x))
  }
  
  bases <- tolower(bases)
  result <- apply(seq.mat, 2, function(b) {
    b <- b[b %in% bases]
    diversity(b)
  })
  names(result) <- 1:length(result)
  result
}