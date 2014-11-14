#' @export summarize.loci
#' 
#' @title Locus Summaries
#' @description Compile standard set of by-locus summaries for diploid data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param analyses character vector specifying analyses to include in summary. 
#'   Can contain "num.alleles", "num.missing", "prop.missing", "allelic.richness", "exptd.het", "obsvd.het", "theta", "hwe.p".
#'   For all available, use "all". Default of \code{NULL} produces all analyses but "hwe.p".
#' @param ... arguments to be passed on to analysis functions.
#' 
#' @return A list with data.frames of locus summary statistics for each stratum with one row per loci. 
#'   The final element in the list, "All", summarizes loci across all individuals in \code{g}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{num.alleles}}, \code{\link{num.missing}},
#'   \code{\link{allelic.richness}}, \code{\link{exptd.het}} \code{\link{obsvd.het}},
#'   \code{\link{theta}}, \code{\link{hwe.genepop}}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#'
#' summarize.loci(msats)

summarize.loci <- function(g, analyses = NULL, ...) {
  stopifnot.gtypes(g, "diploid")
  
  # check analyses
  if(is.null(analyses)) analyses <- c("num.alleles", "num.missing", "prop.missing", "allelic.richness", 
                                      "exptd.het", "obsvd.het", "theta")
  if(!is.character(analyses) | length(analyses) == 0) stop("'analyses' is not character or is zero-length.")
  locus.analyses <- c("num.alleles", "num.missing", "prop.missing", "allelic.richness", 
                      "exptd.het", "obsvd.het", "theta", "hwe.p")
  analyses <- tolower(analyses)
  
  if(analyses[1] == "all") analyses <- locus.analyses
  if(!all(analyses %in% locus.analyses)) stop("some 'analyses' not valid.")
  
  # create analysis order
  analysis.order <- match(analyses, locus.analyses)
  prop.missing <- function(g) num.missing(g, TRUE)
  function.vec <- list(num.alleles, num.missing, prop.missing, 
                       allelic.richness, exptd.het, obsvd.het, theta, hwe.genepop
  )
  
  # run selected functions
  summary.stats <- function(x) {
    result.df <- data.frame(locus = attr(x, "locus.names"))
    function.mat <- sapply(analysis.order, function(i) {
      func <- function.vec[[i]]
      num.args <- length(formals(func))
      if(num.args == 1) func(x) else func(x, ...)
    })
    colnames(function.mat) <- locus.analyses[analysis.order]
    result.df <- cbind(result.df, function.mat)
    rownames(result.df) <- NULL
    result.df
  }
  c(lapply(strata.split(g), summary.stats), list(All = summary.stats(g)))
}