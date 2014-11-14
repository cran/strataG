#' @export stat.phist
#' @importFrom ape dist.dna
#' @importFrom ape as.DNAbin
#' 
#' @title PHIst Analysis of Population Structure
#' 
#' @param g a \code{gtypes} object.
#' @param hap.dist a matrix of pairwise haplotype distances.
#' @param pairwise.deletion passed to \code{\link[ape]{dist.dna}}.
#' @param ... other arguments passed to \code{\link[ape]{dist.dna}} such as 
#'   \code{model}.
#' 
#' @references Excoffier, L., Smouse, P.E., and J.M. Quattro. 1992. Analysis of molecular variance inferred from 
#'   metric distances among DNA haplotypes: Application to human mitochondrial DNA restriction data 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.haps)
#' mtdna <- gtypes(dolph.strata, id.col = 1, strata.col = 2, locus.col = 4, dna.seq = dolph.haps)
#' 
#' stat.phist(mtdna)

stat.phist <- function(g, hap.dist = NULL, pairwise.deletion = TRUE, ...)  {
  g <- check.gtypes(g, "haploid")
  
  stat.name <- "PHIst"
  if(!is.null(hap.dist)) {
    if(!("dist" %in% class(hap.dist))) stop("'hap.dist' must be of class 'dist'")
    if(!all(decode.loci(g) %in% labels(hap.dist))) stop("Some haplotypes in 'g' could not be found in 'hap.dist'")
    if(all(hap.dist[lower.tri(hap.dist)] == 1)) stat.name <- "Fst"
  } else if(is.null(g$sequences)) {
    haps <- unique(g$genotypes[, 2])
    hap.dist <- matrix(1, nrow = length(haps), ncol = length(haps), dimnames = list(haps, haps))
    diag(hap.dist) <- 0
    stat.name <- "Fst"
  } else {
    stopifnot.aligned(g$sequences)
    hap.dist <- dist.dna(as.DNAbin(g$sequences), pairwise.deletion = pairwise.deletion, ...)
  }
  if(!is.matrix(hap.dist)) hap.dist <- as.matrix(hap.dist)
  
# inserted old R code
  
  hap.data <- g$genotypes[g$genotypes[, 1] != -1, , drop = FALSE]
  # Extract summary values
  strata.freq <- table(hap.data[, "strata"])
  num.strata <- length(strata.freq)
  num.samples <- nrow(hap.data)
  strata.hap.freq <- table(hap.data[, "haplotype"], hap.data[, "strata"])
  
  # Calculate sums of squares within strata (Eqn 8a)
  ssd.wp <- sum(sapply(names(strata.freq), function(s) {
    hap.freq <- strata.hap.freq[, s, drop = F]
    hap.freq <- hap.freq[hap.freq[, 1] > 0, , drop = F]
    sum(sapply(rownames(hap.freq), function(h1) {
      freq.1 <- hap.freq[h1, ]
      sapply(rownames(hap.freq), function(h2) hap.dist[h1, h2] * freq.1 * hap.freq[h2, ])
    })) / (2 * strata.freq[s])
  }))
  
  # Calculate sums of squares among strata (Eqn 8b)
  ssd.ap <- sum(sapply(names(strata.freq), function(s1) {
    sapply(names(strata.freq), function(s2) {
      hap.freq.1 <- strata.hap.freq[, s1, drop = F]
      hap.freq.1 <- hap.freq.1[hap.freq.1[, 1] > 0, , drop = F]
      hap.freq.2 <- strata.hap.freq[, s2, drop = F]
      hap.freq.2 <- hap.freq.2[hap.freq.2[, 1] > 0, , drop = F]
      sum(sapply(rownames(hap.freq.1), function(h1) {
        freq.1 <- hap.freq.1[h1, ]
        sapply(rownames(hap.freq.2), function(h2) hap.dist[h1, h2] * freq.1 * hap.freq.2[h2, ])
      }))
    })
  })) / sum(2 * strata.freq)
  ssd.ap <- ssd.ap - ssd.wp
  
  # Calculate average sample size correction for among strata variance 
  #  Eqn 9a in paper, but modified as in Table 8.2.1.1 from Arlequin v3.5.1 manual
  #  (denominator is sum{I} - 1)
  n <- (num.samples - sum(strata.freq ^ 2 / num.samples)) / (num.strata - 1)
  
  # Calculate variance components (Table 1)
  #   Set MSD (SSD / df) equal to expected MSD
  Vc <- ssd.wp / (num.samples - num.strata)
  Vb <- ((ssd.ap / (num.strata - 1)) - Vc) / n
  
  est <- Vb / (Vb + Vc)
  
# lines for C code  
#   est <- if(any(is.na(hap.dist)) | any(is.nan(hap.dist))) NA else {
#     phist_C(g$genotypes[, 2], g$genotypes[, 1], hap.dist)
#   }

  if(is.nan(est)) est <- NA
  
  result <- list(stat.name = stat.name, estimate = est, strata.freq = table(decode.strata(g)))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.phist) <- c(class(stat.phist), "gtype.struct.func")