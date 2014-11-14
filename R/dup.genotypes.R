#' @name dup.genotypes
#' @export dup.genotypes
#' 
#' @title Duplicate Genotypes
#' 
#' @param g a \code{\link[strataG]{gtypes}} object.
#' @param label label for output folder and prefix for files.
#' @param num.shared either number of loci or percentage of loci two individuals must share 
#' to be considered duplicate individuals.
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#' 
#' @return a data.frame with the following columns:
#' \tabular{ll}{
#'   \code{id1, id2} \tab sample ids.\cr
#'   \code{strata1, strata2} \tab sample strata.\cr
#'   \code{num.loci.genotyped} \tab number of loci genotyped for both samples.\cr
#'   \code{num.loci.shared} \tab number of loci shared between both samples.\cr
#'   \code{percent.loci.shared} \tab percent of loci genotyped for both samples that are shared.\cr
#'   \code{mismatch.loci} \tab loci where the two samples do not match.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

dup.genotypes <- function(g, label, num.shared = 0.8, num.cores = 1) {
  stopifnot.gtypes(g)
  
  #if not already, convert num.shared to %
  if(num.shared > 1) num.shared <- num.shared / num.loci(g$genotypes[, -1, drop = FALSE]) 
    
  shared.locs <- pairwise.shared.loci(g, num.cores = num.cores)
  dup.df <- shared.locs[shared.locs[, "prop.same"] >= num.shared, ]
  if(nrow(dup.df) > 0) {
    dup.df$pct.loci.shared <- dup.df$num.same / dup.df$num.not.missing
    strata <- decode.strata(g)
    dup.df$strata1 <- strata[dup.df$id1]
    dup.df$strata2 <- strata[dup.df$id2]
    dup.df$mismatch.loci <- sapply(1:nrow(dup.df), function(i) {
      num.same <- as.matrix(dup.df[i, attr(g, "locus.names")])[1, ]
      loc.names <- names(which(num.same != 4))
      paste(loc.names, collapse = ", ")
    })
    colnames(dup.df)[c(3:5)] <- c("num.loci.shared", "num.loci.genotyped", "percent.loci.shared")
    dup.df <- dup.df[, c("id1", "id2", "strata1", "strata2", "num.loci.genotyped", "num.loci.shared", "percent.loci.shared", "mismatch.loci")]
  } 
  
  if(nrow(dup.df) > 0) {
    sort.order <- order(dup.df$percent.loci.shared, dup.df$num.loci.shared, 
      rev(dup.df$id1), rev(dup.df$id2), decreasing = TRUE
    )
    dup.df <- dup.df[sort.order, ]
  }
  rownames(dup.df) <- NULL
  
  write.csv(dup.df, paste(label, ".duplicate.samples.csv", sep = ""), row.names = FALSE)

  invisible(dup.df)
}
