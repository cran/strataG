#' @export stat.d.jost
#' @importFrom swfscMisc harmonic.mean
#' 
#' @title Jost's D Estimate of Population Divergence
#' 
#' @param g a \code{gtypes} object.
#' 
#' @details Calculate Dest as presented in Jost 2008, Eqn 13.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.d.jost(msats)

stat.d.jost <- function(g) {
  g <- check.gtypes(g, "diploid")
  
  g <- decode(g)
  terms <- sapply(1:length(attr(g, "locus.names")), function(i) {
    locus <- na.omit(cbind(strata = g$genotypes[, "strata"], g$genotypes[, locus.cols(i)]))
    if(nrow(locus) == 0) return(c(a = NA, b = NA))
    
    allele.mat <- locus[, 1:2]
    allele.mat <- rbind(allele.mat, locus[, c(1, 3)])
    num.strata <- length(unique(locus[, 1]))
    allele.strata.freq <- table(allele.mat[, 2], allele.mat[, 1])
    i.terms <- sapply(rownames(allele.strata.freq), function(allele) {
      j.terms <- sapply(colnames(allele.strata.freq), function(strata) {
        Nj <- sum(allele.strata.freq[, strata])
        Nij <- allele.strata.freq[allele, strata]
        a.term1 <- Nij / Nj
        a.term2 <- a.term1 ^ 2
        b.term <- Nij * (Nij - 1) / (Nj * (Nj - 1))
        c(a.term1 = a.term1, a.term2 = a.term2, b.term = b.term)
      })
      a.term1 <- sum(j.terms["a.term1", ]) ^ 2
      a.term2 <- sum(j.terms["a.term2", ])
      c(a = (a.term1 - a.term2) / (num.strata - 1), b = sum(j.terms["b.term", ]))
    })
    
    c(a = sum(i.terms["a", ]), b = sum(i.terms["b", ]))
  })
  est <- 1 - terms["a", ] / terms["b", ]
  
  result <- list(stat.name = "D (Jost 2008)", estimate = harmonic.mean(est), strata.freq = table(g$genotypes[, "strata"]))
  class(result) <- c(class(result), "gtype.struct.stat")
  result
}
class(stat.d.jost) <- c(class(stat.d.jost), "gtype.struct.func")