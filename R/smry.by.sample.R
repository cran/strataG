#' @name smry.by.sample
#' @export smry.by.sample
#' 
#' @title By-Sample Summaries
#' 
#' @param g a \code{\link[strataG]{gtypes}} object.
#' @param label label for output folder and prefix for files.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

smry.by.sample <- function(g, label) {
  stopifnot.gtypes(g)
  g <- decode(g)
  
  number.loci <- num.loci(g$genotypes[, -1, drop = FALSE])
  by.sample.smry <- do.call(rbind, lapply(1:nrow(g$genotypes), function(id) {
    sample <- rownames(g$genotypes)[id]
    
    num.missing <- sum(sapply(1:number.loci, function(i) {
      locus <- rbind(g$genotypes[id, as.vector(locus.cols(i))])
      sum(apply(locus, 1, function(x) any(is.na(x))))
    }))
    
    homozygous <- sapply(1:number.loci, function(i) {
      locus <- na.omit(g$genotypes[id, locus.cols(i), drop = FALSE])
      if(length(locus) == 0) return(NA)
      ifelse(locus[, 1] == locus[, 2], 1, 0)
    })
    
    c(sample = sample, 
      num.loci.missing.genotypes = num.missing,
      pct.loci.missing.genotypes = num.missing / number.loci,
      pct.loci.homozygous = mean(homozygous, na.rm = TRUE)
    )
  }))                 
  write.csv(by.sample.smry, paste(label, ".sample.summary.csv", sep = ""), row.names = FALSE)
  
  invisible(NULL)
}
