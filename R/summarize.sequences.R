#' @export summarize.sequences
#' 
#' @title DNA Sequence Summaries
#' @description Generate a standard set of summaries for sequences.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @note The function labels haplotypes first so may run slower for large data sets with individual sequences.

#' @return A data.frame with sequence summary statistics with one row per stratum and a final row
#'   for the all sequences in \code{g}. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{haplotypic.diversity}}, \code{\link{variable.sites}}

summarize.sequences <- function(g) {
  stopifnot.gtypes(g, "haploid")
  g <- label.haplotypes(g)
  
  summary.stats <- function(x) {
    num.samples <- nrow(x$genotypes)
    haps <- as.character(unique(x$genotypes[, 2]))
    num.haps <- length(haps)
    het <- haplotypic.diversity(x)
    num.Ns <- sapply(haps, function(h) sum(tolower(x$sequences[[h]])  == "n"))
    vs <- variable.sites(x$sequences, c("a", "c", "g", "t", "-")) 
    c(num.samples = num.samples, num.haps = num.haps, num.var.sites = ncol(vs$site.freqs),
      het = het, median.num.Ns = median(num.Ns), max.num.Ns = max(num.Ns)
    )
  }
  
  t(cbind(sapply(strata.split(g, remove.sequences = TRUE), summary.stats), All = summary.stats(g)))
}
