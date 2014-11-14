#' @name smry.by.locus
#' @export smry.by.locus
#' 
#' @title By-Locus Summaries
#' 
#' @param g a \code{\link[strataG]{gtypes}} object.
#' @param label label for output folder and prefix for files.
#' 
#' @note Requires that GENEPOP is installed on the system and accessible from the command line.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

smry.by.locus <- function(g, label) {
  stopifnot.gtypes(g)
  
  by.locus.smry <- function(g, file) {
    loci.num <- num.loci(g$genotypes[, -1, drop = FALSE])
    result.summary <- do.call(rbind, lapply(attr(g, "locus.names"), function(loc) {
      subset.g <- subset(g, loci = loc) 
      num.genotyped <- nrow(g$genotypes) - num.missing(subset.g)
      c(loc, nrow(g$genotypes), num.genotyped, num.genotyped / nrow(g$genotypes), num.alleles(subset.g),
        allelic.richness(subset.g), exptd.het(subset.g), obsvd.het(subset.g), theta(subset.g),
        hwe.genepop(subset.g)
      )
    }))
    colnames(result.summary) <- c("locus", "num.samples", "num.genotyped", "pct.genotyped", 
     "num.alleles", "allelic.richness", "expected.heterozygosity", "observed.heterozygosity", "theta", "hwe"
    )
    write.csv(result.summary, file = file, row.names = FALSE)      
  }
  
  # Summaries for each stratum
  strata.g <- strata.split(g)
  for(x in names(strata.g)) {
    file <- paste(label, "locus.summary", x, "csv", sep = ".")
    by.locus.smry(strata.g[[x]], file) 
  }

  # Summaries for all strata
  file <- paste(label, ".locus.summary.all.strata.csv", sep = "")
  by.locus.smry(g, file)
  
  invisible(NULL)
}
