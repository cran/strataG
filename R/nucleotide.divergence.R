#' @export nucleotide.divergence
#' @importFrom ape dist.dna
#' @importFrom ape as.DNAbin
#' 
#' @title Nucleotide Divergence
#' @description Calculate Nei's dA between strata, and distributions of 
#'   between- and within-strata nucleotide divergence (sequence distance).
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param probs a numeric vector of probabilities of the pairwise distance distributions with 
#'   values in \code{[0,1]}.
#' @param ... arguments passed to \code{\link[ape]{dist.dna}} such as 
#'   \code{model} or \code{pairwise.deletion}.
#' 
#' @return a list with summaries of the \code{within} and \code{between} strata 
#'   pairwise distances including Nei's dA. 
#'   
#' @references Nei, M., and S. Kumar (2000) Molecular Evolution and Phylogenetics. 
#'   Oxford University Press, Oxford. (dA: pp. 256, eqn 12.67)
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

nucleotide.divergence <- function(g, probs = c(0, 0.025, 0.5, 0.975, 1), ...) {  
  stopifnot.gtypes(g, "haploid")
  stopifnot.aligned(g$sequences)  
  
  pair.dist.summary <- function(sample.pairs, hap.dist) {
    pairwise.dist <- apply(sample.pairs, 1, function(haps) {
      if(any(is.na(haps))) NA else as.vector(hap.dist[haps[1], haps[2]])
    })
    dist.quant <- quantile(pairwise.dist, probs, na.rm = TRUE)
    names(dist.quant) <- paste("pct.", probs, sep = "")
    c(mean = mean(pairwise.dist, na.rm = TRUE), dist.quant)
  }
  
  gd <- decode(g)  
  hap.dist <- dist.dna(as.DNAbin(gd$sequences), as.matrix = T, ...)
  strata <- sort(unique(gd$genotypes[, 1]))
  within.dist <- t(sapply(strata, function(x) {
    gd.haps <- gd$genotypes[gd$genotypes[, 1] == x, 2]
    pair.dist.summary(t(combn(gd.haps, 2)), hap.dist)
  }))
  
  strata.pairs <- t(combn(strata, 2))
  between.dist <- t(apply(strata.pairs, 1, function(sp) {
    haps.1 <- gd$genotypes[gd$genotypes[, 1] == sp[1], 2]
    haps.2 <- gd$genotypes[gd$genotypes[, 1] == sp[2], 2]
    result <- pair.dist.summary(expand.grid(haps.1, haps.2), hap.dist)
    dA <- result["mean"] - (sum(within.dist[sp, "mean"], na.rm = TRUE) / 2)
    names(dA) <- NULL
    c(dA = dA, result)
  }))
  between.dist <- data.frame(strata.pairs, between.dist)
  colnames(between.dist)[1:2] <- c("strata.1", "strata.2")  
  
  list(within = within.dist, between = between.dist)
}