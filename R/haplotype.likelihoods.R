#' @export haplotype.likelihoods
#' @importFrom MASS fitdistr
#' @importFrom ape dist.dna
#' @importFrom ape as.DNAbin
#' 
#' @title Haplotype Likelihoods
#' @description Calculate likelihood of each haplotype based on gamma distribution of pairwise distances.
#' 
#' @param dna.seq a list of DNA sequences.
#' @param model a character string specifying the evolutionary model to be used. Passed to 
#'   \code{\link[ape]{dist.dna}}. 
#' @param pairwise.deletion a logical indicating whether to delete the sites with missing data 
#'   in a pairwise way. Passed to \code{\link[ape]{dist.dna}}. 
#' 
#' @return vector of delta(log-Likelihoods) for each haplotype, sorted from smallest to largest, and 
#'   a plot of their distributions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(dolph.haps)
#' haplotype.likelihoods(dolph.haps)

haplotype.likelihoods <- function(dna.seq, model = "raw", pairwise.deletion = FALSE) {
  stopifnot.aligned(dna.seq)
  
  seq.dist <- dist.dna(as.DNAbin(dna.seq), as.matrix = TRUE, model = model, 
    pairwise.deletion = pairwise.deletion
  )
  if(any(is.nan(seq.dist)) | any(is.na(seq.dist))) {
    warning("NA/NaN in pairwise distance matrix. NULL returned.")
    return(NULL)
  }
  dist.vec <- seq.dist[lower.tri(seq.dist)]
  gamma.coeffs <- fitdistr(dist.vec, "gamma", lower = 0)$estimate
  
  log.lik <- sapply(rownames(seq.dist), function(x) {
    x.dist <- seq.dist[x, ]
    x.dist <- x.dist[names(x.dist) != x]
    sum(log(dgamma(x.dist, gamma.coeffs[1], gamma.coeffs[2])), na.rm = TRUE)
  })
  delta.log.lik <- sort(log.lik - max(log.lik, na.rm = T), decreasing = FALSE)
  
  dotchart(rev(delta.log.lik), pch = 19, bg = "black", 
    xlab = expression(paste(Delta, "lnL")), main = "Haplotype Likelihoods"
  )
  
  delta.log.lik
}