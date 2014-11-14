#' @export nei.Da
#' @aliases Da
#' 
#' @title Nei's Da
#' @description Calcuate frequency-based Nei's Da for haploid or diploid data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param na.rm logical. Delete loci which have missing data for every sample in a stratum?
#' 
#' @details Returns Nei's Da for each pair of strata.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Nei et al 1983 Accuracy of Estimated Phylogenetic Trees from 
#'   Molecular Data. J Mol Evol 19:153-170 (eqn 7)\cr
#'   Nei, M., and S. Kumar (2000) Molecular Evolution and Phylogenetics. 
#'   Oxford University Press, Oxford. (pp. 268, eqn 13.6)

nei.Da <- function(g, na.rm = FALSE) {
  stopifnot.gtypes(g)
  
  strata <- unique(decode.strata(g))
  if(length(strata) == 1) stop("More than one stratum required to calculate Nei's Da.") 
  strata.pairs <- t(combn(sort(strata), 2))

  Da <- if(is.diploid(g)) {
    apply(strata.pairs, 1, function(s) {
      pair.gtypes <- subset.gtypes(g, strata = s)
      loc.sum <- sapply(attr(pair.gtypes, "locus.names"), function(x) {
        locus <- as.vector(pair.gtypes$genotypes[, locus.cols(x, pair.gtypes)])
        locus[locus == -1] <- NA
        freqs <- prop.table(table(locus, rep(pair.gtypes$genotypes[, "strata"], 2)))
        sum(apply(freqs, 1, function(f) if(all(f == 0)) NA else sqrt(prod(f))))
      })
      1 - sum(loc.sum, na.rm = na.rm) / sum(!is.na(loc.sum))
    })
  } else {
    apply(strata.pairs, 1, function(s) {
      pair.gtypes <- subset.gtypes(g, strata = s)
      freqs <- prop.table(hap.freqs(pair.gtypes))
      1 - sum(apply(freqs, 1, function(f) if(all(f == 0)) NA else sqrt(prod(f))), na.rm = na.rm)
    })  
  }
  result <- as.data.frame(strata.pairs)
  result <- cbind(result, Da)
  colnames(result) <- c("strata.1", "strata.2", "Da")
  result
}