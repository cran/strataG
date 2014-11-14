#' @export allele.freqs
#' 
#' @title Allele Frequencies
#' 
#' @param g a \code{gtypes} object.
#' @param use.na.rows logical. FALSE deletes any row with an NA for either allele. 
#'   TRUE uses all unique (non-NA) alleles in frequency.
#' 
#' @return A list of allele or haplotype frequencies.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#'
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' 
#' allele.freqs(msats)[c(1, 5)]

allele.freqs <- function(g, use.na.rows = FALSE) {  
  stopifnot.gtypes(g)
  
  if(is.haploid(g)) return(hap.freqs(g))

  g <- decode(g)
  sapply(attr(g, "locus.names"), function(x) {
    i <- which(attr(g, "locus.names") == x)
    locus <- g$genotypes[, locus.cols(i)]
    allele <- if(use.na.rows) as.vector(na.omit(locus)) else na.omit(as.vector(locus))
    allele.tbl <- table(allele)
    result <- data.frame(locus = rep(x, length(allele.tbl)))
    result$allele <- names(allele.tbl)
    result$freq <- as.vector(allele.tbl)
    result$prop <- result$freq / sum(result$freq)
    result
  }, USE.NAMES = TRUE, simplify = FALSE)
}