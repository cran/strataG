#' @rdname pop.diff.test
#' @export pop.diff.test pairwise.test overall.test
#' 
#' @title Population Differentiation Tests
#' @description Conduct overall and/or pairwise tests of population differentiation.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param stats a character vector specifying which anlayses to conduct. For haploid \code{g}, can contain 
#'   "fst", "phist", or "chi2", or for diploid \code{g}, can contain "fst", "fst.prime", "gst", "gst.prime", 
#'   "gst.dbl.prime", "d", or "chi2". The default is "all" which executes all available tests for given ploidy.
#' @param overall logical. Calculate overall statistic of population differentiation?
#' @param pairwise logical. Calculate statistics of population differentiation between all pairs of strata?
#' @param nrep number specifying number of permutation replicates to use for permutation test.
#' @param keep.null logical. Keep the null distribution from the permutation test?
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#' @param quietly logical. Print progress output to screen?
#' @param write.output logical. Write a .csv file with results?
#' @param ... other parameters to be passed to population differentiation functions.
#' 
#' @return list with results
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' # Conduct an overall test of differentiation
#' overall.test(msats, stats = "chi2", nrep = 20, write.output = FALSE)
#' 
#' # Conduct a pairwise test between all strata
#' pairwise.test(msats, stats = c("gst.prime", "gst.dbl.prime"), nrep = 20, write.output = FALSE)
#' 
#' # Conduct both overall and pairwise tests
#' result <- pop.diff.test(msats, stats = c("fst", "fst.prime"), nrep = 20, write.output = FALSE)
#' print(result$overall)
#' print(result$pairwise)

pop.diff.test <- function(g, stats = "all", overall = TRUE, pairwise = TRUE, nrep = 1000, 
                           keep.null = FALSE, num.cores = 1, quietly = FALSE, write.output = FALSE, ...) {
  two.strata <- length(unique(decode.strata(g))) == 2
  
  overall.results <- if(overall | two.strata) {
    overall.test(g, stats = stats, nrep = nrep, keep.null = keep.null, 
      num.cores = num.cores, quietly = quietly, write.output = write.output, ...
    )
  } else NULL
  
  pairwise.results <- if(pairwise & !two.strata) {
    pairwise.test(g, stats = stats, nrep = nrep, keep.null = keep.null, 
      num.cores = num.cores, quietly = quietly, write.output = write.output, ...
    )
  } else NULL
  
  invisible(list(overall = overall.results, pairwise = pairwise.results))
}